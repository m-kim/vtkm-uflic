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
#include <fstream>
#include <string>

#include "Evaluator.h"
#include "Integrator.h"
#include "UFLIC.h"
#include "Draw.h"


#include <iostream>
#include <fstream>
#include <chrono>

typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

void saveAs(std::string fileName, 
	vtkm::cont::ArrayHandle<vtkm::Float32 > canvasArray, vtkm::Id Width, vtkm::Id Height) {
	std::ofstream of(fileName.c_str(), std::ios_base::binary | std::ios_base::out);
	of << "P6" << std::endl << Width << " " << Height << std::endl << 255 << std::endl;
	//ColorBufferType::PortalConstControl colorPortal = this->ColorBuffer.GetPortalConstControl();
	for (vtkm::Id yIndex = Height - 1; yIndex >= 0; yIndex--)
	{
		for (vtkm::Id xIndex = 0; xIndex < Width; xIndex++)
		{
			vtkm::Float32 val = canvasArray.GetPortalConstControl().Get(yIndex * Width + xIndex);
			val *= 255;
			vtkm::Vec<vtkm::Float32, 4> tuple(val, val, val, val);
			of << (unsigned char)(tuple[0]);
			of << (unsigned char)(tuple[1]);
			of << (unsigned char)(tuple[2]);
		}
	}
	of.close();
}


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
  IntegratorType integrator(eval, 2.0);


  ParticleAdvectionWorkletType advect(integrator);
  advect.Run(sl, sr, fieldArray);


  //for (vtkm::Id i = 0; i < sr.GetNumberOfValues(); i++)
  //{
  //  vtkm::Vec<FieldType, Size> p = sr.GetPortalConstControl().Get(i);
  //  std::cout << p[0] << " " << p[1] << std::endl;
  //}


  vtkm::cont::DataSetBuilderUniform dataSetBuilder;
  vtkm::cont::DataSet ds = dataSetBuilder.Create(dim);

	std::vector<FieldType> canvas[2], omega(dim[0] * dim[1], 0);
	canvas[0].resize(dim[0] * dim[1]);
	canvas[1].resize(dim[1] * dim[0]);

	for (int i = 0; i < canvas[0].size(); i++) {
		canvas[0][i] =  rand() / (float)RAND_MAX;
		canvas[1][i] = 0;

	}
  vtkm::cont::ArrayHandle<FieldType > canvasArray[2], omegaArray;
  canvasArray[0] = vtkm::cont::make_ArrayHandle(&canvas[0][0], canvas[0].size());
	canvasArray[1] = vtkm::cont::make_ArrayHandle(&canvas[1][0], canvas[1].size());
	omegaArray = vtkm::cont::make_ArrayHandle(&omega[0], omega.size());
  DrawLineWorkletType drawline(ds);
  drawline.Run(canvasArray[0], canvasArray[1], omegaArray,sl,sr);


	saveAs("uflic.pnm", canvasArray[1], dim[0], dim[1]);
	//vtkm::rendering::Mapper mapper;
	//vtkm::rendering::Canvas canvas(512, 512);
	//vtkm::rendering::Scene scene;

	//scene.AddActor(vtkm::rendering::Actor(
	//	ds.GetCellSet(), ds.GetCoordinateSystem(), ds.GetField(fieldNm), colorTable));
	//vtkm::rendering::Camera camera;
	//SetCamera<ViewType>(camera, ds.GetCoordinateSystem().GetBounds());

	//vtkm::rendering::View2D view;
}
