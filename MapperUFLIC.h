#ifndef MAPPERUFLIC_H
#define MAPPERUFLIC_H
#include <vtkm/rendering/Mapper.h>
#include <vtkm/rendering/raytracing/RayTracer.h>
#include <vtkm/rendering/raytracing/Logger.h>
#include <vtkm/rendering/raytracing/TriangleExtractor.h>
#include <vtkm/rendering/raytracing/RayOperations.h>
#include <vtkm/rendering/raytracing/Ray.h>
#include <vtkm/cont/Timer.h>
#include "CanvasUFLIC.h"
namespace detail
{

class RayStatusFilter : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  RayStatusFilter() {}
  using ControlSignature = void(FieldIn<>, FieldInOut<>);
  using ExecutionSignature = void(_1, _2);
  VTKM_EXEC
  void operator()(const vtkm::Id& hitIndex, vtkm::UInt8& rayStatus) const
  {
    if (hitIndex == -1)
      rayStatus = RAY_EXITED_DOMAIN;
    else if (rayStatus != RAY_EXITED_DOMAIN && rayStatus != RAY_TERMINATED)
      rayStatus = RAY_ACTIVE;
    //else printf("Bad status state %d \n",(int)rayStatus);
  }
}; //class RayStatusFileter

class RayMapCanvas : public vtkm::worklet::WorkletMapField
{
protected:
  vtkm::Matrix<vtkm::Float32, 4, 4> InverseProjView;
  vtkm::Id Width;
  vtkm::Float32 DoubleInvHeight;
  vtkm::Float32 DoubleInvWidth;
  vtkm::Vec<vtkm::Float32, 3> Origin;

public:
  VTKM_CONT
  RayMapCanvas(const vtkm::Matrix<vtkm::Float32, 4, 4>& inverseProjView,
               const vtkm::Id width,
               const vtkm::Id height,
               const vtkm::Vec<vtkm::Float32, 3>& origin)
    : InverseProjView(inverseProjView)
    , Width(width)
    , Origin(origin)
  {
    VTKM_ASSERT(width > 0);
    VTKM_ASSERT(height > 0);
    DoubleInvHeight = 2.f / static_cast<vtkm::Float32>(height);
    DoubleInvWidth = 2.f / static_cast<vtkm::Float32>(width);
  }

  using ControlSignature = void(FieldIn<>, FieldInOut<>, WholeArrayIn<>);
  using ExecutionSignature = void(_1, _2, _3);

  template <typename Precision, typename DepthPortalType>
  VTKM_EXEC void operator()(const vtkm::Id& pixelId,
                            Precision& maxDistance,
                            const DepthPortalType& depths) const
  {
    vtkm::Vec<vtkm::Float32, 4> position;
    position[0] = static_cast<vtkm::Float32>(pixelId % Width);
    position[1] = static_cast<vtkm::Float32>(pixelId / Width);
    position[2] = static_cast<vtkm::Float32>(depths.Get(pixelId));
    position[3] = 1;
    // transform into normalized device coordinates (-1,1)
    position[0] = position[0] * DoubleInvWidth - 1.f;
    position[1] = position[1] * DoubleInvHeight - 1.f;
    position[2] = 2.f * position[2] - 1.f;
    // offset so we don't go all the way to the same point
    position[2] -= 0.00001f;
    position = vtkm::MatrixMultiply(InverseProjView, position);
    vtkm::Vec<vtkm::Float32, 3> p;
    p[0] = position[0] / position[3];
    p[1] = position[1] / position[3];
    p[2] = position[2] / position[3];
    p = p - Origin;

    maxDistance = vtkm::Magnitude(p);
  }
}; //class RayMapMinDistances

} // namespace detail
void MapCanvasToRays(vtkm::rendering::raytracing::Ray<vtkm::Float32>& rays,
                                    const vtkm::rendering::Camera& camera,
                                    const vtkm::rendering::Canvas &canvas)
{
  vtkm::Id width = canvas.GetWidth();
  vtkm::Id height = canvas.GetHeight();
  vtkm::Matrix<vtkm::Float32, 4, 4> projview =
    vtkm::MatrixMultiply(camera.CreateProjectionMatrix(width, height), camera.CreateViewMatrix());
  bool valid;
  vtkm::Matrix<vtkm::Float32, 4, 4> inverse = vtkm::MatrixInverse(projview, valid);
  if (!valid)
    throw vtkm::cont::ErrorBadValue("Inverse Invalid");

  vtkm::worklet::DispatcherMapField<detail::RayMapCanvas>(
    detail::RayMapCanvas(inverse, width, height, camera.GetPosition()))
    .Invoke(rays.PixelIdx, rays.MaxDistance, canvas.GetDepthBuffer());
}

class VTKM_RENDERING_EXPORT MapperUFLIC : public vtkm::rendering::Mapper
{
public:
    struct InternalsType
    {
      CanvasUFLIC* Canvas;
      vtkm::rendering::raytracing::RayTracer Tracer;
      vtkm::rendering::raytracing::Camera RayCamera;
      vtkm::rendering::raytracing::Ray<vtkm::Float32> Rays;
      bool CompositeBackground;
      bool Shade;
      VTKM_CONT
      InternalsType()
        : Canvas(nullptr)
        , CompositeBackground(true)
        , Shade(true)
      {
      }
    };
    MapperUFLIC()
      : Internals(new InternalsType)
    {
    }

    ~MapperUFLIC()
    {
    }

    void SetCanvas(vtkm::rendering::Canvas* canvas) override
    {
      if (canvas != nullptr)
      {
        this->Internals->Canvas = dynamic_cast<CanvasUFLIC*>(canvas);
        if (this->Internals->Canvas == nullptr)
        {
          throw vtkm::cont::ErrorBadValue("Ray Tracer: bad canvas type. Must be CanvasRayTracer");
        }
      }
      else
      {
        this->Internals->Canvas = nullptr;
      }
    }

    vtkm::rendering::Canvas* GetCanvas() const override
    {
      return this->Internals->Canvas;
    }

    void RenderCells(const vtkm::cont::DynamicCellSet& cellset,
                                      const vtkm::cont::CoordinateSystem& coords,
                                      const vtkm::cont::Field& scalarField,
                                      const vtkm::cont::ColorTable& vtkmNotUsed(colorTable),
                                      const vtkm::rendering::Camera& camera,
                                      const vtkm::Range& scalarRange) override
    {
      vtkm::rendering::raytracing::Logger* logger = vtkm::rendering::raytracing::Logger::GetInstance();
      logger->OpenLogEntry("mapper_ray_tracer");
      vtkm::cont::Timer<> tot_timer;
      vtkm::cont::Timer<> timer;

      // make sure we start fresh
      this->Internals->Tracer.Clear();
      //
      // Add supported shapes
      //
      vtkm::Bounds shapeBounds;
      vtkm::rendering::raytracing::TriangleExtractor triExtractor;
      triExtractor.ExtractCells(cellset);
      if (triExtractor.GetNumberOfTriangles() > 0)
      {
        vtkm::rendering::raytracing::TriangleIntersector* triIntersector = new vtkm::rendering::raytracing::TriangleIntersector();
        triIntersector->SetData(coords, triExtractor.GetTriangles());
        this->Internals->Tracer.AddShapeIntersector(triIntersector);
        shapeBounds.Include(triIntersector->GetShapeBounds());
      }

      //
      // Create rays
      //
      vtkm::rendering::raytracing::Camera& cam = this->Internals->Tracer.GetCamera();
      cam.SetParameters(camera, *this->Internals->Canvas);
      this->Internals->RayCamera.SetParameters(camera, *this->Internals->Canvas);

      this->Internals->RayCamera.CreateRays(this->Internals->Rays, shapeBounds);
      this->Internals->Rays.Buffers.at(0).InitConst(0.f);
      MapCanvasToRays(
        this->Internals->Rays, camera, *this->Internals->Canvas);



      this->Internals->Tracer.SetField(scalarField, scalarRange);

      this->Internals->Tracer.SetColorMap(this->ColorMap);
      this->Internals->Tracer.SetShadingOn(this->Internals->Shade);
      this->Internals->Tracer.Render(this->Internals->Rays);

      timer.Reset();
      this->Internals->Canvas->WriteToCanvas(
        this->Internals->Rays, this->Internals->Rays.Buffers.at(0).Buffer, camera);

      if (this->Internals->CompositeBackground)
      {
        this->Internals->Canvas->BlendBackground();
      }

      vtkm::Float64 time = timer.GetElapsedTime();
      logger->AddLogData("write_to_canvas", time);
      time = tot_timer.GetElapsedTime();
      logger->CloseLogEntry(time);
    }

    void SetCompositeBackground(bool on)
    {
      this->Internals->CompositeBackground = on;
    }

    void SetShadingOn(bool on)
    {
      this->Internals->Shade = on;
    }

    void StartScene() override
    {
      // Nothing needs to be done.
    }

    void EndScene() override
    {
      // Nothing needs to be done.
    }

    vtkm::rendering::Mapper* NewCopy() const override
    {
      return new MapperUFLIC(*this);
    }

private:
  std::shared_ptr<InternalsType> Internals;

  struct RenderFunctor;

};

#endif
