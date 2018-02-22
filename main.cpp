#include <vtkm/Math.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/rendering/View2D.h>
#include <vtkm/rendering/Canvas.h>
#include <vtkm/rendering/Mapper.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/exec/RandomArray.h>
#include <fstream>
#include <sstream>
#include <memory>

#include "Evaluator.h"
#include "Integrator.h"
#include "UFLIC.h"
#include "Draw.h"
#include "Sharpen.h"
#include "Normalize.h"
#include "Jitter.h"
#include "Convolve.h"

#include "Reader.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <type_traits>


bool do_print = false;

typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

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

class SetRandomArray : public vtkm::worklet::WorkletMapField
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

template<typename VecComponentType>
void saveAs(std::string fileName,
  vtkm::cont::ArrayHandle<VecComponentType > canvasArray,
  vtkm::Id Width, vtkm::Id Height) {
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


template<typename EvalType, typename VecType, vtkm::IdComponent Size>
void run( std::shared_ptr<Reader<VecType, Size>> reader)
{
  typedef vtkm::Int32 FieldType;
  typedef vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>> VecHandle;
  typedef RK4Integrator<EvalType, VecType, Size> IntegratorType;
  typedef ParticleAdvectionWorklet<IntegratorType, VecType, Size, DeviceAdapter> ParticleAdvectionWorkletType;
  typedef DrawLineWorklet<FieldType, VecType, Size, DeviceAdapter> DrawLineWorkletType;

  const vtkm::IdComponent ttl = 4;

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
  vtkm::worklet::DispatcherMapField<ResetParticles<VecType, Size>, DeviceAdapter> resetDispatcher(dim[0]);
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


  vtkm::cont::ArrayHandle<FieldType > canvasArray[ttl], propFieldArray[2], omegaArray, texArray;
  for (int i = 0; i < ttl; i++) {
    canvasArray[i].Allocate(dim[0] * dim[1]);
  }
  propFieldArray[0].Allocate(dim[0] * dim[1]);
  propFieldArray[1].Allocate(dim[0] * dim[1]);

  omegaArray.Allocate(dim[0] * dim[1]);
  texArray.Allocate(dim[0] * dim[1]);

  randomDispatcher.Invoke(indexArray, texArray);
  DrawLineWorkletType drawline(bounds, dim);
  DoNormalize<FieldType, DeviceAdapter> donorm(dim);
  DoSharpen<FieldType, DeviceAdapter> dosharp(dim);
  DoJitter<FieldType, DeviceAdapter> dojitter(dim);

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

    vtkm::worklet::DispatcherMapField<zero_voxel, DeviceAdapter> zeroDispatcher;
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

    t += dt;// / (vtkm::Float32)ttl + 1.0 / (vtkm::Float32)ttl;

  }

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

std::tuple<int,int,int>
parse(int argc, char **argv){

  int x = 512;
  int y = 256;
  int which = 0;

  for (int i=1; i<argc; i++){
    if (!strcmp(argv[i], "bfield")){
      which = 1;
    }
    else if (!strcmp(argv[i], "PSI")){
      which = 2;
    }
    else if (!strcmp(argv[i], "dims")){
      if (i+1 < argc && i+2 < argc){
        x = atoi(argv[i+1]);
        y = atoi(argv[i+2]);
        i += 2;
      }
    }
    else if (!strcmp(argv[i], "print")){
      do_print = true;
    }
  }

  return std::make_tuple(which,x,y);
}

int main(int argc, char **argv)
{
  const int Size = 2;
  typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;
  typedef vtkm::Float32 VecType;

  std::tuple<int, int, int> ret;

  ret = parse(argc,argv);

  std::shared_ptr<Reader<VecType, Size>> reader;

  if (std::get<0>(ret) == 1){
    reader = std::shared_ptr<ReaderVTK<VecType, Size>>(new ReaderVTK<VecType, Size>("/home/mkim/vtkm-uflic/BField_2d.vtk", 12));
    run<VectorField<VecType,Size>,VecType,Size>(reader);
  }

  else if (std::get<0>(ret) == 2){
    //std::shared_ptr<Reader<VecType, Size,  ReaderPS<VecType, Size,ReaderXGC<VecType,Size>>>> reader(new ReaderPS<VecType, Size, ReaderXGC<VecType,Size>>("/home/mkim/vtkm-uflic/psi2q/2D_packed/psi2D_packed_normalized_256_99.vec", vtkm::Id2(256,256), Bounds(0,256,0,256)));
    reader = std::shared_ptr<ReaderXGC<VecType, Size>>(new ReaderXGC<VecType, Size>("/home/mkim/vtkm-uflic/psi2q/2D_packed/psi2D_packed_512_", vtkm::Id2(512,512), Bounds(0,512,0,512), 12));
    run<VectorField<VecType,Size>,VecType,Size>(reader);
  }
  else{

    int x = std::get<1>(ret);
    int y = std::get<2>(ret);
    reader = std::shared_ptr<ReaderCalc<VecType, Size>>(new ReaderCalc<VecType, Size>("XGC_", vtkm::Id2(x,y), Bounds(0,x,0,y), vtkm::Vec<VecType,Size>(2,1), 12));
    run<DoubleGyreField<VecType,Size>,VecType,Size>(reader);
  }
}
