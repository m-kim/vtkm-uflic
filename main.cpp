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
#include "Sharpen.h"


#include <iostream>
#include <fstream>
#include <chrono>

typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

template<typename VecComponentType>
void saveAs(std::string fileName, 
	vtkm::cont::ArrayHandle<VecComponentType > canvasArray, 
	vtkm::cont::ArrayHandle<VecComponentType> omegaArray,
	vtkm::Id Width, vtkm::Id Height) {
	std::ofstream of(fileName.c_str(), std::ios_base::binary | std::ios_base::out);
	of << "P6" << std::endl << Width << " " << Height << std::endl << 255 << std::endl;
	//ColorBufferType::PortalConstControl colorPortal = this->ColorBuffer.GetPortalConstControl();
	for (vtkm::Id yIndex = Height - 1; yIndex >= 0; yIndex--)
	{
		for (vtkm::Id xIndex = 0; xIndex < Width; xIndex++)
		{
			VecComponentType val = canvasArray.GetPortalConstControl().Get(yIndex * Width + xIndex);
			VecComponentType omega = omegaArray.GetPortalConstControl().Get(yIndex * Width + xIndex);
			if (omega > 0)
				val /= omega;
			vtkm::Vec<VecComponentType, 4> tuple(val, val, val, val);
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
  typedef vtkm::Float32 VecType;
	typedef vtkm::Int32 FieldType;

  typedef vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>> VecHandle;
  typedef VecHandle::template ExecutionTypes<DeviceAdapter>::PortalConst VecPortalConstType;

  typedef DoubleGyreField<VecPortalConstType, VecType> EvalType;
  typedef EulerIntegrator<EvalType, VecType, Size> IntegratorType;

  typedef ParticleAdvectionWorklet<IntegratorType, VecType, Size, DeviceAdapter> ParticleAdvectionWorkletType;
  typedef DrawLineWorklet<FieldType, VecType, Size, DeviceAdapter> DrawLineWorkletType;


  vtkm::Id2 dim(256,256);

  std::vector<vtkm::Vec<VecType, Size>> pl, pr;

  for (int y=0;y<dim[0]; y++){
    for (int x=0; x<dim[1]; x++){
      pl.push_back(vtkm::Vec<VecType,Size>(x+0.5, y+0.5));
      pr.push_back(vtkm::Vec<VecType,Size>(x+0.5, y+0.5));
    }
  }
  vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>> sl, sr;
  sl = vtkm::cont::make_ArrayHandle(pl);
  sr = vtkm::cont::make_ArrayHandle(pr);

  vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>> VecArray;
  EvalType eval(Bounds(0,dim[1],0,dim[1]));
  IntegratorType integrator(eval, 3.0);
	

  ParticleAdvectionWorkletType advect(integrator);
  advect.Run(sl, sr, VecArray);


  //for (vtkm::Id i = 0; i < sr.GetNumberOfValues(); i++)
  //{
  //  vtkm::Vec<VecType, Size> p = sr.GetPortalConstControl().Get(i);
  //  std::cout << p[0] << " " << p[1] << std::endl;
  //}


  vtkm::cont::DataSetBuilderUniform dataSetBuilder;
  vtkm::cont::DataSet ds = dataSetBuilder.Create(dim);

	std::vector<FieldType> canvas[2], omega(dim[0] * dim[1], 0);
	canvas[0].resize(dim[0] * dim[1]);
	canvas[1].resize(dim[1] * dim[0]);

	for (int i = 0; i < canvas[0].size(); i++) {
		canvas[0][i] = rand() % RAND_MAX;
		canvas[1][i] = 0;
	}

  vtkm::cont::ArrayHandle<FieldType > canvasArray[2], omegaArray;
  canvasArray[0] = vtkm::cont::make_ArrayHandle(&canvas[0][0], canvas[0].size());
	canvasArray[1] = vtkm::cont::make_ArrayHandle(&canvas[1][0], canvas[1].size());
	omegaArray = vtkm::cont::make_ArrayHandle(&omega[0], omega.size());
  DrawLineWorkletType drawline(ds);
  drawline.Run(canvasArray[0], canvasArray[1], omegaArray,sl,sr);

	DoSharpen<FieldType, DeviceAdapter> dosharp(dim);
	dosharp.Run(canvasArray[1], canvasArray[0]);
	


	saveAs("uflic.pnm", canvasArray[0], omegaArray, dim[0], dim[1]);
	//vtkm::rendering::Mapper mapper;
	//vtkm::rendering::Canvas canvas(512, 512);
	//vtkm::rendering::Scene scene;

	//scene.AddActor(vtkm::rendering::Actor(
	//	ds.GetCellSet(), ds.GetCoordinateSystem(), ds.GetField(fieldNm), colorTable));
	//vtkm::rendering::Camera camera;
	//SetCamera<ViewType>(camera, ds.GetCoordinateSystem().GetBounds());

	//vtkm::rendering::View2D view;
}
