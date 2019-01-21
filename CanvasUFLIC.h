#ifndef CANVAS_UFLIC_H
#define CANVAS_UFLIC_H
#include <vtkm/rendering/Canvas.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/rendering/raytracing/Ray.h>

namespace myinternal
{
class MySurfaceConverter : public vtkm::worklet::WorkletMapField
{
  vtkm::Matrix<vtkm::Float32, 4, 4> ViewProjMat;

public:
  VTKM_CONT
  MySurfaceConverter(const vtkm::Matrix<vtkm::Float32, 4, 4> viewProjMat)
    : ViewProjMat(viewProjMat)
  {
  }

  using ControlSignature =
    void(FieldIn<>, WholeArrayInOut<>, FieldIn<>, FieldIn<>, FieldIn<>, WholeArrayOut<>, WholeArrayOut<>);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, WorkIndex);
  template <typename Precision,
            typename ColorPortalType,
            typename DepthBufferPortalType,
            typename ColorBufferPortalType>
  VTKM_EXEC void operator()(const vtkm::Id& pixelIndex,
                            ColorPortalType& colorBufferIn,
                            const Precision& inDepth,
                            const vtkm::Vec<Precision, 3>& origin,
                            const vtkm::Vec<Precision, 3>& dir,
                            DepthBufferPortalType& depthBuffer,
                            ColorBufferPortalType& colorBuffer,
                            const vtkm::Id& index) const
  {
    vtkm::Vec<Precision, 3> intersection = origin + inDepth * dir;
    vtkm::Vec<vtkm::Float32, 4> point;
    point[0] = static_cast<vtkm::Float32>(intersection[0]);
    point[1] = static_cast<vtkm::Float32>(intersection[1]);
    point[2] = static_cast<vtkm::Float32>(intersection[2]);
    point[3] = 1.f;

    vtkm::Vec<vtkm::Float32, 4> newpoint;
    newpoint = vtkm::MatrixMultiply(this->ViewProjMat, point);
    newpoint[0] = newpoint[0] / newpoint[3];
    newpoint[1] = newpoint[1] / newpoint[3];
    newpoint[2] = newpoint[2] / newpoint[3];

    vtkm::Float32 depth = newpoint[2];

    if (depth > depthBuffer.Get(index))
        return;
    depth = 0.5f * (depth) + 0.5f;
    vtkm::Vec<vtkm::Float32, 4> color;
    color[0] = static_cast<vtkm::Float32>(colorBufferIn.Get(index * 4 + 0));
    color[1] = static_cast<vtkm::Float32>(colorBufferIn.Get(index * 4 + 1));
    color[2] = static_cast<vtkm::Float32>(colorBufferIn.Get(index * 4 + 2));
    color[3] = static_cast<vtkm::Float32>(colorBufferIn.Get(index * 4 + 3));

    vtkm::Vec<vtkm::Float32, 4> vec_point = newpoint;
    vec_point += color;
    vec_point[3] = 1.f;

    vtkm::Vec<vtkm::Float32, 4> new_vec_point;
    new_vec_point = vtkm::MatrixMultiply(this->ViewProjMat, vec_point);
    vtkm::Vec<vtkm::Float32,3> fin_vec_point;
    new_vec_point[0] = new_vec_point[0]/ new_vec_point[3];
    new_vec_point[1] = new_vec_point[1]/ new_vec_point[3];
    new_vec_point[2] = new_vec_point[2]/ new_vec_point[3];
    new_vec_point = new_vec_point - newpoint;
    new_vec_point[3] = 1.f;

    fin_vec_point[0] = new_vec_point[0];
    fin_vec_point[1] = new_vec_point[1];
    fin_vec_point[2] = new_vec_point[2];
    fin_vec_point = vtkm::Normal(fin_vec_point);

    new_vec_point[0] = fin_vec_point[0];
    new_vec_point[1] = fin_vec_point[1];
    new_vec_point[2] = fin_vec_point[2];
    new_vec_point[3] = 1.0f;
//     // blend the mapped color with existing canvas color
//     vtkm::Vec<vtkm::Float32, 4> inColor = colorBuffer.Get(pixelIndex);

//     // if transparency exists, all alphas have been pre-multiplied
//     vtkm::Float32 alpha = (1.f - color[3]);
//     color[0] = color[0] + inColor[0] * alpha;
//     color[1] = color[1] + inColor[1] * alpha;
//     color[2] = color[2] + inColor[2] * alpha;
//     color[3] = inColor[3] * alpha + color[3];

//     // clamp
//     for (vtkm::Int32 i = 0; i < 4; ++i)
//     {
//       color[i] = vtkm::Min(1.f, vtkm::Max(color[i], 0.f));
//     }
//     // The existing depth should already been feed into the ray mapper
//     // so no color contribution will exist past the existing depth.

    depthBuffer.Set(pixelIndex, depth);
    colorBuffer.Set(pixelIndex, new_vec_point);
  }
}; //class SurfaceConverter
template <typename Precision>
VTKM_CONT void WriteToCanvas(const vtkm::rendering::raytracing::Ray<Precision>& rays,
                             const vtkm::cont::ArrayHandle<Precision>& colors,
                             const vtkm::rendering::Camera& camera,
                             vtkm::Id width,
                             vtkm::Id height,
                             vtkm::rendering::Canvas::DepthBufferType &depthBuffer,
                            vtkm::rendering::Canvas::ColorBufferType &colorBuffer)
{
    std::cout <<"My write to canvas!" << std::endl;
  vtkm::Matrix<vtkm::Float32, 4, 4> viewProjMat =
    vtkm::MatrixMultiply(camera.CreateProjectionMatrix(width, height),
                         camera.CreateViewMatrix());

  vtkm::worklet::DispatcherMapField<MySurfaceConverter>(MySurfaceConverter(viewProjMat))
    .Invoke(rays.PixelIdx,
            colors,
            rays.Distance,
            rays.Origin,
            rays.Dir,
            depthBuffer,
            colorBuffer);

  //Force the transfer so the vectors contain data from device
  colorBuffer.GetPortalControl().Get(0);
  depthBuffer.GetPortalControl().Get(0);
}
} // namespace internal

class VTKM_RENDERING_EXPORT CanvasUFLIC : public vtkm::rendering::Canvas
{
public:
  CanvasUFLIC(vtkm::Id width = 1024, vtkm::Id height = 1024)
  : vtkm::rendering::Canvas(width,height){std::cout <<"CanvasUFLIC" << std::endl;}

  ~CanvasUFLIC(){}

  vtkm::rendering::Canvas* NewCopy() const
{
  return new CanvasUFLIC(*this);
}

  void WriteToCanvas(const vtkm::rendering::raytracing::Ray<vtkm::Float32>& rays,
                                    const vtkm::cont::ArrayHandle<vtkm::Float32>& colors,
                                    const vtkm::rendering::Camera& camera)
{
std::cout << "writetocanvas 1" << std::endl;
  myinternal::WriteToCanvas(rays, colors, camera, this->GetWidth(), this->GetHeight(),
                            this->GetDepthBuffer(), this->GetColorBuffer());
}

void WriteToCanvas(const vtkm::rendering::raytracing::Ray<vtkm::Float64>& rays,
                                    const vtkm::cont::ArrayHandle<vtkm::Float64>& colors,
                                    const vtkm::rendering::Camera& camera)
{
std::cout << "writetocanvas 2" << std::endl;
  myinternal::WriteToCanvas(rays, colors, camera, this->GetWidth(), this->GetHeight(),
                            this->GetDepthBuffer(), this->GetColorBuffer());
}

}; // class CanvasRayTracer

#endif