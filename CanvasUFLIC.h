#ifndef CANVAS_UFLIC_H
#define CANVAS_UFLIC_H
#include <vtkm/rendering/Canvas.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/rendering/raytracing/Ray.h>
#include <vtkm/cont/ArrayCopy.h>

namespace myinternal
{
class MySurfaceConverter : public vtkm::worklet::WorkletMapField
{
  vtkm::Matrix<vtkm::Float32, 4, 4> ViewProjMat;
  vtkm::Id width, height;
public:
  VTKM_CONT
  MySurfaceConverter(const vtkm::Matrix<vtkm::Float32, 4, 4> viewProjMat,
                    vtkm::Id w,
                    vtkm::Id h)
    : ViewProjMat(viewProjMat),
    width(w),
    height(h)
  {
  }

  using ControlSignature =
    void(FieldIn<>, FieldOut<>, WholeArrayInOut<>, FieldIn<>, FieldIn<>, FieldIn<>, WholeArrayOut<>, WholeArrayOut<>);
  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8, WorkIndex);
  template <typename Precision,
            typename ColorPortalType,
            typename DepthBufferPortalType,
            typename ColorBufferPortalType>
  VTKM_EXEC void operator()(const vtkm::Id& pixelIndex,
                            vtkm::Vec<vtkm::Float32,2> &vel,
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
    
    depth = 0.5f * (depth) + 0.5f;
    vtkm::Vec<vtkm::Float32, 4> color;
    color[0] = static_cast<vtkm::Float32>(colorBufferIn.Get(index * 4 + 0));
    color[1] = static_cast<vtkm::Float32>(colorBufferIn.Get(index * 4 + 1));
    color[2] = static_cast<vtkm::Float32>(colorBufferIn.Get(index * 4 + 2));
    color[3] = static_cast<vtkm::Float32>(colorBufferIn.Get(index * 4 + 3));
    
    if (vtkm::Magnitude(color) > 1e-6){

      vtkm::Vec<vtkm::Float32, 4> inColor = colorBuffer.Get(pixelIndex);

      vtkm::Vec<vtkm::Float32, 4> newPos = color;
      newPos[0] = intersection[0] + newPos[0];
      newPos[1] = intersection[1] + newPos[1];
      newPos[2] = intersection[2] + newPos[2];
      newPos[3] = 1.0;


      auto newcolor = vtkm::MatrixMultiply(this->ViewProjMat, newPos);
      newcolor[0] = newcolor[0] / newcolor[3];
      newcolor[1] = newcolor[1] / newcolor[3];
      newcolor[2] = newcolor[2] / newcolor[3];

      newcolor[0] = 0.5f * newcolor[0] + 0.5f;
      newcolor[1] = 0.5f * newcolor[1] + 0.5f;
      newcolor[2] = 0.5f * newcolor[2] + 0.5f;
      newcolor[0] *= width;
      newcolor[1] *= height;
      newcolor[0] = vtkm::Max(0, newcolor[0]);
      newcolor[0] = vtkm::Min(width, newcolor[0]);
      newcolor[1] = vtkm::Max(0, newcolor[1]);
      newcolor[1] = vtkm::Min(height, newcolor[1]);


      vtkm::Vec<vtkm::Float32, 2> pxpos(pixelIndex % width, pixelIndex / width);
      pxpos[0] = vtkm::Max(0, pxpos[0]);
      pxpos[0] = vtkm::Min(width, pxpos[0]);
      pxpos[1] = vtkm::Max(0, pxpos[1]);
      pxpos[1] = vtkm::Min(height, pxpos[1]);

      vel[0] = newcolor[0] - pxpos[0];
      vel[1] = newcolor[1] - pxpos[1];
      vtkm::Normalize(vel);

  //    // blend the mapped color with existing canvas color
  //    // if transparency exists, all alphas have been pre-multiplied
  //    vtkm::Float32 alpha = (1.f - color[3]);
  //    color[0] = color[0] + inColor[0] * alpha;
  //    color[1] = color[1] + inColor[1] * alpha;
  //    color[2] = color[2] + inColor[2] * alpha;
  //    color[3] = inColor[3] * alpha + color[3];

      // // clamp
      // for (vtkm::Int32 i = 0; i < 4; ++i)
      // {
      //   color[i] = vtkm::Min(1.f, vtkm::Max(color[i], 0.f));
      // }
      // // The existing depth should already been feed into the ray mapper
      // // so no color contribution will exist past the existing depth.

      depthBuffer.Set(pixelIndex, depth);
      colorBuffer.Set(pixelIndex, color);
    }
    else{
      colorBuffer.Set(pixelIndex, vtkm::Vec<vtkm::Float32,4>(0,0,0,0));
      vel[0] = 0;
      vel[1] = 0;
      depthBuffer.Set(pixelIndex, -1.0);

    }
  }
}; //class SurfaceConverter
template <typename Precision>
VTKM_CONT void WriteToCanvas(const vtkm::rendering::raytracing::Ray<Precision>& rays,
                             const vtkm::cont::ArrayHandle<Precision>& colors,
                             const vtkm::rendering::Camera& camera,
                             vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,2>> &pixelPos,
                             vtkm::Id width,
                             vtkm::Id height,
                             vtkm::rendering::Canvas::DepthBufferType &depthBuffer,
                            vtkm::rendering::Canvas::ColorBufferType &colorBuffer)
{
  vtkm::Matrix<vtkm::Float32, 4, 4> viewProjMat =
    vtkm::MatrixMultiply(camera.CreateProjectionMatrix(width, height),
                         camera.CreateViewMatrix());

  vtkm::worklet::DispatcherMapField<MySurfaceConverter>(MySurfaceConverter(viewProjMat, width,height))
    .Invoke(rays.PixelIdx,
            pixelPos,
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
  using DepthBufferType = vtkm::cont::ArrayHandle<vtkm::Float32>;

  CanvasUFLIC(vtkm::Id width = 1024, vtkm::Id height = 1024)
  : vtkm::rendering::Canvas(width,height),
    iter_cnt(1)
  {
  }

  ~CanvasUFLIC(){}

  vtkm::rendering::Canvas* NewCopy() const
{
  return new CanvasUFLIC(*this);
}

  void WriteToCanvas(const vtkm::rendering::raytracing::Ray<vtkm::Float32>& rays,
                                    const vtkm::cont::ArrayHandle<vtkm::Float32>& colors,
                                    const vtkm::rendering::Camera& camera,
                                    const vtkm::Float32 stepsize)
{
  pixelPos.Allocate(rays.PixelIdx.GetNumberOfValues());
  vtkm::cont::ArrayHandleConstant<vtkm::Id> zero(0, pixelPos.GetNumberOfValues());
  vtkm::cont::ArrayCopy(zero, pixelPos);

  myinternal::WriteToCanvas(rays, colors, camera, pixelPos, this->GetWidth(), this->GetHeight(),
                            this->GetDepthBuffer(), this->GetColorBuffer());

  pixelVel.Allocate(this->GetWidth() * this->GetHeight());
  vtkm::cont::ArrayCopy(
              vtkm::cont::ArrayHandleConstant<vtkm::Vec<vtkm::Float32,2>>(vtkm::Vec<vtkm::Float32,2>(0), pixelVel.GetNumberOfValues()),
              pixelVel);

  pixelIdx.Allocate(pixelPos.GetNumberOfValues());
  vtkm::cont::ArrayCopy(rays.PixelIdx, pixelIdx);

  for (int i=0; i < rays.PixelIdx.GetNumberOfValues(); i++){
    vtkm::Id id = rays.PixelIdx.GetPortalConstControl().Get(i);

    auto vel = pixelPos.GetPortalControl().Get(i);
    pixelVel.GetPortalControl().Set(id, vel);
  }

}

//void WriteToCanvas(const vtkm::rendering::raytracing::Ray<vtkm::Float64>& rays,
//                                    const vtkm::cont::ArrayHandle<vtkm::Float64>& colors,
//                                    const vtkm::rendering::Camera& camera)
//{
//  pixelPos.Allocate(rays.PixelIdx.GetNumberOfValues());
//  vtkm::cont::ArrayHandleConstant<vtkm::Vec<vtkm::Int32,2>> zero(0, pixelPos.GetNumberOfValues());
//  vtkm::cont::ArrayCopy(zero, pixelPos);

//  myinternal::WriteToCanvas(rays, colors, camera, pixelPos, this->GetWidth(), this->GetHeight(),
//                            this->GetDepthBuffer(), this->GetColorBuffer());
//}

vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,2>> pixelPos, pixelVel;
vtkm::cont::ArrayHandle<vtkm::Id> pixelIdx;
int iter_cnt;

}; // class CanvasRayTracer

#endif
