#include <vtkm/Math.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include "Evaluator.h"
#include "Integrator.h"
#include "UFLIC.h"
#include "Draw.h"


#include <iostream>
#include <fstream>
#include <chrono>

typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

//std::shared_ptr<vtkm::cont::DataSetBuilderExplicitIterative> dataSetBuilder;


int main(int argc, char **argv)
{


  const int Size = 2;
  typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;
  typedef vtkm::Float32 FieldType;

  typedef vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>> FieldHandle;
  typedef FieldHandle::template ExecutionTypes<DeviceAdapter>::PortalConst FieldPortalConstType;

  typedef DoubleGyreField<FieldPortalConstType, FieldType> EvalType;
  typedef EulerIntegrator<EvalType, FieldType, Size> IntegratorType;

  typedef ParticleAdvectionWorklet<IntegratorType, FieldType, Size, DeviceAdapter> ParticleAdvectionWorkletType;
  typedef DrawLineWorklet<FieldType, Size, DeviceAdapter> DrawLineWorkletType;


  vtkm::Id2 dim(256,256);

  std::vector<vtkm::Vec<FieldType, Size>> pl, pr;

  for (int y=0;y<dim[0]; y++){
    for (int x=0; x<dim[1]; x++){
      pl.push_back(vtkm::Vec<FieldType,Size>(x+0.5, y+0.5));
      pr.push_back(vtkm::Vec<FieldType,Size>(x+0.5, y+0.5));
    }
  }
  vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>> sl, sr;
  sl = vtkm::cont::make_ArrayHandle(pl);
  sr = vtkm::cont::make_ArrayHandle(pr);

  vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>> fieldArray;
  EvalType eval(Bounds(0,dim[1],0,dim[1]));
  IntegratorType integrator(eval, 0.1);


  ParticleAdvectionWorkletType advect(integrator);
  advect.Run(sl, sr, fieldArray);


  for (vtkm::Id i = 0; i < sr.GetNumberOfValues(); i++)
  {
    vtkm::Vec<FieldType, Size> p = sr.GetPortalConstControl().Get(i);
    std::cout << p[0] << " " << p[1] << std::endl;
  }


  vtkm::cont::DataSetBuilderUniform dataSetBuilder;
  vtkm::cont::DataSet ds = dataSetBuilder.Create(dim);

  std::vector<FieldType> canvas(dim[0]*dim[1],0), omega(dim[0]*dim[1],0), noise(dim[0]*dim[1],0);

  vtkm::cont::ArrayHandle<FieldType > canvasArray, omegaArray, noiseArray;
  canvasArray = vtkm::cont::make_ArrayHandle(&canvas[0], canvas.size());
  omegaArray = vtkm::cont::make_ArrayHandle(&omega[0], omega.size());
  noiseArray = vtkm::cont::make_ArrayHandle(&noise[0], noise.size());
  DrawLineWorkletType drawline(ds);
  drawline.Run(canvasArray, omegaArray,sl,sr);

}
