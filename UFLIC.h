#ifndef UFLIC_H
#define UFLIC_H
#include <chrono>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/BinaryOperators.h>
#include "Evaluator.h"
#include "Integrator.h"
#include "ParticleAdvection.h"
#include "Draw.h"
#include "Sharpen.h"
#include "Normalize.h"
#include "Jitter.h"
#include "Convolve.h"
#include "RandomMapField.h"
#include "RandomArray.h"
#include "Reader.h"
#include "Quiver.h"

template<typename FieldType, vtkm::Id Size>
class ResetParticles : public vtkm::worklet::WorkletMapField
{
public:
  typedef  vtkm::Vec<FieldType, Size> VecType;
  typedef void ControlSignature(FieldIn<>, FieldOut<>);
  typedef void ExecutionSignature(_1, _2);
  //
  VTKM_CONT
  ResetParticles(vtkm::Id d,
                 vtkm::Id2 s = vtkm::Vec<vtkm::Id,2>(1,1))
    : dim(d)
    , span(s) {}

  VTKM_EXEC_CONT void operator()(const vtkm::Id&idx,
                                 VecType& out) const
  {
    vtkm::Id y = (idx / dim) * span[1];
    vtkm::Id x = (idx % dim) * span[0];
    out = VecType(x, y);
  }

  vtkm::Id dim;
  vtkm::Id2 span;
};

class SetRandomArray : public RandomWorklet
{
public:
  typedef void ControlSignature(FieldIn<>, RandomArrayInOut<>);
  typedef void ExecutionSignature(WorkIndex, _2);
  typedef _1 InputDomain;

  VTKM_CONT
  SetRandomArray(vtkm::Vec<vtkm::Int32,2> ab)
    :a(ab[0]),
      b(ab[1])
  {

  }
  template <typename RandomArrayType>
  VTKM_EXEC void operator()(const vtkm::Id& index, const RandomArrayType& randomArray) const
  {
    typedef typename RandomArrayType::ValueType ValueType;
    randomArray.Set(static_cast<ValueType>(index),a,b);
  }

private:
  vtkm::Int32 a,b;
};

template<typename EvalType, typename VecType, vtkm::IdComponent Size>
class UFLIC
{
public:
    typedef vtkm::Int32 FieldType;
    typedef vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>> VecHandle;
    typedef RK4Integrator<EvalType, VecType, Size> IntegratorType;
    typedef ParticleAdvectionWorklet<IntegratorType, VecType, Size> ParticleAdvectionWorkletType;
    typedef DrawLineWorklet<FieldType, VecType, Size> DrawLineWorkletType;
    using ArrayType = vtkm::cont::ArrayHandle<FieldType >;
    UFLIC(int _ttl = 4)
    : ttl(_ttl)
    , do_print(false)
    {
      canvasArray.resize(ttl);
    }
    void run( std::shared_ptr<Reader<VecType, Size>> reader)
    {
      const vtkm::Id loop_cnt = reader->iter_cnt;
      reader->readFile();

      #ifdef VTKM_CUDA
      cudaFree(0);
      #endif
      auto t0 = std::chrono::high_resolution_clock::now();

      vtkm::Id2 dim = reader->dim;
      vtkm::Vec<VecType,Size> spacing = reader->spacing;
      Bounds bounds = reader->bounds;

      vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>> vecArray;


      std::vector<vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>>> sl(ttl), sr(ttl);
      vtkm::worklet::DispatcherMapField<ResetParticles<VecType, Size>> resetDispatcher(dim[0]);
      vtkm::worklet::DispatcherMapField<SetRandomArray> randomDispatcher(vtkm::Vec<vtkm::Int32, 2>(0, 255));

      vtkm::cont::ArrayHandleConstant<vtkm::Int8> mask(1, dim[0]*dim[1]);

      vtkm::cont::ArrayHandleCounting<vtkm::Id> indexArray(vtkm::Id(0), 1, dim[0]*dim[1]);

      for (int i = 0; i < ttl; i++) {
          sl[i].Allocate(dim[0] * dim[1]);
          resetDispatcher.Invoke(indexArray, sl[i]);
          sr[i].Allocate(dim[0] * dim[1]);
          resetDispatcher.Invoke(indexArray, sr[i]);

      }

      vtkm::Float32 t = 0;
      const vtkm::Float32 dt = 0.1;


      for (int i = 0; i < ttl; i++) {
          canvasArray[i].Allocate(dim[0] * dim[1]);
      }
      propFieldArray[0].Allocate(dim[0] * dim[1]);
      propFieldArray[1].Allocate(dim[0] * dim[1]);

      omegaArray.Allocate(dim[0] * dim[1]);
      texArray.Allocate(dim[0] * dim[1]);

      randomDispatcher.Invoke(indexArray, texArray);
      DrawLineWorkletType drawline(dim);
      DoNormalize<FieldType> donorm(dim);
      DoSharpen<FieldType> dosharp(dim);
      DoJitter<FieldType> dojitter(dim);

      for (int loop = 0; loop < loop_cnt; loop++) {
          EvalType eval(t, Bounds(0, dim[0], 0, dim[1]), spacing);
          IntegratorType integrator(eval, 3.0);
          ParticleAdvectionWorkletType advect(integrator);

          reader->next(vecArray);
          resetDispatcher.Invoke(indexArray, sl[loop%ttl]);
          //reset the current canvas
      #ifdef VTKM_CUDA
          randomDispatcher.Invoke(indexArray, canvasArray[loop%ttl]);
      #else
          for (int i=0; i<canvasArray[loop%ttl].GetNumberOfValues(); i++){
          canvasArray[loop%ttl].GetPortalControl().Set(i,rand()%255);
          }
      #endif
          vtkm::cont::ArrayHandleConstant<FieldType> zero(0, propFieldArray[0].GetNumberOfValues());
          vtkm::cont::ArrayCopy(zero, propFieldArray[0]);
          vtkm::cont::ArrayCopy(zero, propFieldArray[1]);
          vtkm::cont::ArrayCopy(zero, omegaArray);

          for (int i = 0; i < vtkm::Min(ttl, loop+1); i++) {
          advect.Run(sl[i], sr[i], vecArray);
          drawline.Run(canvasArray[i], propFieldArray[0], omegaArray, mask, indexArray, sl[i], sr[i]);
          }

          sr.swap(sl);


          donorm.Run(propFieldArray[0], omegaArray, propFieldArray[1]);
          if (do_print){
          std::stringstream fn;
          fn << "uflic-" << loop << ".pnm";
          saveAs(fn.str().c_str(), propFieldArray[1], dim[0], dim[1]);
          }

          //REUSE omegaArray as a temporary cache to sharpen
          dosharp.Run(propFieldArray[1], omegaArray);
          dojitter.Run(omegaArray, texArray, canvasArray[(loop) % ttl]);

          //t += dt;// / (vtkm::Float32)ttl + 1.0 / (vtkm::Float32)ttl;

      }

      result = propFieldArray[1];
      auto t1 = std::chrono::high_resolution_clock::now();

      std::cout << "Finished dt: " << dt << " cnt: " << loop_cnt << " time: " << std::chrono::duration<double>(t1-t0).count() << " s" << std::endl;
      std::stringstream fn;
      fn << "uflic-final" << ".pnm";
      saveAs(fn.str().c_str(), propFieldArray[1], dim[0], dim[1]);


      //vtkm::rendering::Mapper mapper;
      //vtkm::rendering::Canvas canvas(512, 512);
      //vtkm::rendering::Scene scene;

      //scene.AddActor(vtkm::rendering::Actor(
      //	ds.GetCellSet(), ds.GetCoordinateSystem(), ds.GetField(fieldNm), colorTable));
      //vtkm::rendering::Camera camera;
      //SetCamera<ViewType>(camera, ds.GetCoordinateSystem().GetBounds());

      //vtkm::rendering::View2D view;

    }
    template<typename VecComponentType>
    static void saveAs(std::string fileName,
                vtkm::cont::ArrayHandle<VecComponentType > canvasArray,
                vtkm::Id Width, vtkm::Id Height) 
    {
      std::ofstream of(fileName.c_str(), std::ios_base::binary | std::ios_base::out);
      of << "P6" << std::endl << Width << " " << Height << std::endl << 255 << std::endl;
      //ColorBufferType::PortalConstControl colorPortal = this->ColorBuffer.GetPortalConstControl();
      for (vtkm::Id yIndex = Height - 1; yIndex >= 0; yIndex--)
      {
          for (vtkm::Id xIndex = 0; xIndex < Width; xIndex++)
          {
          VecComponentType val = canvasArray.GetPortalConstControl().Get(yIndex * Width + xIndex);

          vtkm::Vec<VecComponentType, 4> tuple(val, val, val, val);
          of << (unsigned char)(tuple[0]);
          of << (unsigned char)(tuple[1]);
          of << (unsigned char)(tuple[2]);
          }
      }
      of.close();
    }

    template<typename VecArrayType>
    static void saveQuiver(
        vtkm::Id2 dim,
        VecArrayType &vec
        )
    {
      vtkm::Id2 spacing(8,8);
      vtkm::cont::ArrayHandleCounting<vtkm::Id> indexArray(vtkm::Id(0), vtkm::Id(1), dim[0]/spacing[0]*dim[1]/spacing[1]);

      vtkm::worklet::DispatcherMapField<ResetParticles<vtkm::Float32, 2>> resetDispatcher(ResetParticles<vtkm::Float32,2>(dim[0]/spacing[0], spacing));
      typename std::remove_reference<decltype(vec)>::type pl;
      pl.Allocate(dim[0]/spacing[0]*dim[1]/spacing[1]);
      resetDispatcher.Invoke(indexArray, pl);


      saveQuiver(dim, pl, vec, spacing);
    }
    template<typename VecArrayType>
    static void saveQuiver(
        vtkm::Id2 dim,
        VecArrayType &pl,
        VecArrayType &vec,
        vtkm::Id2 spacing
        )
    {
      using CanvasArrayType = typename vtkm::cont::ArrayHandle<FieldType>;
      Quiver<VecArrayType, CanvasArrayType, vtkm::Id2> q(dim);
      CanvasArrayType canvas;
      canvas.Allocate(dim[0] * dim[1]);
      q.draw(pl, vec, canvas, spacing);
      saveAs("output.pnm", canvas, dim[0], dim[1]);

    }
    bool do_print;
    std::vector<ArrayType> canvasArray;
    ArrayType propFieldArray[2],result, omegaArray, texArray;
    const vtkm::IdComponent ttl;

};
#endif
