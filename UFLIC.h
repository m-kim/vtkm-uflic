#ifndef UFLIC_H
#define UFLIC_H
#include <chrono>
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

class zero_voxel : public vtkm::worklet::WorkletMapField
{
public:
  typedef void ControlSignature(FieldIn<>, FieldOut<>);
  typedef void ExecutionSignature(_1, WorkIndex, _2);
  //
  VTKM_CONT
  zero_voxel() {}

  template <typename T>
  VTKM_EXEC_CONT void operator()(const vtkm::Id&,
                                 const vtkm::Id& vtkmNotUsed(index),
                                 T& voxel_value) const
  {
    voxel_value = T(0);
  }
};

template<typename FieldType, vtkm::Id Size>
class ResetParticles : public vtkm::worklet::WorkletMapField
{
public:
  typedef  vtkm::Vec<FieldType, Size> VecType;
  typedef void ControlSignature(FieldIn<>, FieldOut<>);
  typedef void ExecutionSignature(_1, _2);
  //
  VTKM_CONT
  ResetParticles(vtkm::Id &d) : dim(d) {}

  VTKM_EXEC_CONT void operator()(const vtkm::Id&idx,
                                 VecType& out) const
  {
    vtkm::Id y = idx / dim;
    vtkm::Id x = idx % dim;
    out = VecType(x + 0.5, y + 0.5);
  }

  vtkm::Id dim;
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
      DrawLineWorkletType drawline(bounds, dim);
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

          vtkm::worklet::DispatcherMapField<zero_voxel> zeroDispatcher;
          zeroDispatcher.Invoke(indexArray, propFieldArray[0]);
          zeroDispatcher.Invoke(indexArray, propFieldArray[1]);
          zeroDispatcher.Invoke(indexArray, omegaArray);

          for (int i = 0; i < vtkm::Min(ttl, loop+1); i++) {
          advect.Run(sl[i], sr[i], vecArray);
          drawline.Run(canvasArray[i], propFieldArray[0], omegaArray, sl[i], sr[i]);
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
    void saveAs(std::string fileName,
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

    bool do_print;
    std::vector<ArrayType> canvasArray;
    ArrayType propFieldArray[2],result, omegaArray, texArray;
    const vtkm::IdComponent ttl;

};
#endif
