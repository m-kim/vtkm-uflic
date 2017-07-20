#include <vtkm/Math.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include "Evaluator.h"
#include "Integrator.h"
#include "UFLIC.h"

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



  std::vector<vtkm::Vec<FieldType, Size>> pl, pr;

  for (int y=0;y<256; y++){
    for (int x=0; x<256; x++){
      pl.push_back(vtkm::Vec<FieldType,Size>(x+0.5, y+0.5));
      pr.push_back(vtkm::Vec<FieldType,Size>(x+0.5, y+0.5));
    }
  }
  vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>> sl, sr;
  sl = vtkm::cont::make_ArrayHandle(pl);
  sr = vtkm::cont::make_ArrayHandle(pr);

  vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>> fieldArray;
  EvalType eval(vtkm::Bounds(0,256,0,256,0,0));
  IntegratorType integrator(eval, 0.1);


  ParticleAdvectionWorkletType worklet(integrator);
  worklet.Run(sl, sr, fieldArray);

  for (vtkm::Id i = 0; i < sr.GetNumberOfValues(); i++)
  {
    vtkm::Vec<FieldType, Size> p = sr.GetPortalConstControl().Get(i);
    std::cout << p[0] << " " << p[1] << std::endl;
  }

}
