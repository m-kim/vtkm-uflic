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
#include "Normalize.h"

#include <iostream>
#include <fstream>
#include <chrono>

typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

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


int main(int argc, char **argv)
{


  const int Size = 2;
  typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;
  typedef vtkm::Float32 VecType;
	typedef vtkm::Int32 FieldType;

  typedef vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>> VecHandle;
  typedef VecHandle::template ExecutionTypes<DeviceAdapter>::PortalConst VecPortalConstType;

  typedef DoubleGyreField<VecPortalConstType, VecType> EvalType;
  typedef RK4Integrator<EvalType, VecType, Size> IntegratorType;

  typedef ParticleAdvectionWorklet<IntegratorType, VecType, Size, DeviceAdapter> ParticleAdvectionWorkletType;
  typedef DrawLineWorklet<FieldType, VecType, Size, DeviceAdapter> DrawLineWorkletType;


  const vtkm::Id2 dim(256,256);
	const vtkm::IdComponent slice = 4;
  std::vector<vtkm::Vec<VecType, Size>> pl[slice], pr[slice];

	for (int i = 0; i < slice; i++) {
		for (int y = 0; y<dim[0]; y++) {
			for (int x = 0; x<dim[1]; x++) {
				pl[i].push_back(vtkm::Vec<VecType, Size>(x + 0.5, y + 0.5));
				pr[i].push_back(vtkm::Vec<VecType, Size>(x + 0.5, y + 0.5));
			}
		}
	}
	vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>> sl[slice], sr[slice];
	for (int i = 0; i<slice; i++) {
		sl[i] = vtkm::cont::make_ArrayHandle(pl[i]);
		sr[i] = vtkm::cont::make_ArrayHandle(pr[i]);
	}


	vtkm::cont::DataSetBuilderUniform dataSetBuilder;
	vtkm::cont::DataSet ds = dataSetBuilder.Create(dim);

	std::vector<FieldType> canvas[slice], propertyField[2], omega(dim[0] * dim[1], 0);
	vtkm::Float32 t = 0;
	const vtkm::Float32 dt = 0.1;
	
	for (int i = 0; i < 2; i++) {
		propertyField[i].resize(dim[0] * dim[1], 0);
	}

	for (int i = 0; i < slice; i++) {
		canvas[i].resize(dim[0] * dim[1], 0);
	}

	for (int i = 0; i < canvas[0].size(); i++) {
		canvas[0][i] = rand() % 256;
	}

	vtkm::cont::ArrayHandle<FieldType > canvasArray[slice], propFieldArray[2], omegaArray;
	for (int i = 0; i < slice; i++) {
		canvasArray[i] = vtkm::cont::make_ArrayHandle(&canvas[i][0], canvas[i].size());
	}
	propFieldArray[0] = vtkm::cont::make_ArrayHandle(&propertyField[0][0], propertyField[0].size());
	propFieldArray[1] = vtkm::cont::make_ArrayHandle(&propertyField[1][0], propertyField[1].size());
	omegaArray = vtkm::cont::make_ArrayHandle(&omega[0], omega.size());

  vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>> VecArray;

	for (int loop = 0; loop < 1; loop++) {
		for (int i = 0; i < propFieldArray[0].GetNumberOfValues(); i++) {
			propFieldArray[0].GetPortalControl().Set(i, 0);
			propFieldArray[1].GetPortalControl().Set(i, 0);
			omegaArray.GetPortalControl().Set(i, 0);
		}
		for (int i = 0; i < 1; i++) {
			EvalType eval(t, Bounds(0, dim[0], 0, dim[1]));
			IntegratorType integrator(eval, 3.0);
			ParticleAdvectionWorkletType advect(integrator);
			DrawLineWorkletType drawline(ds);


			advect.Run(sl[i], sr[i], VecArray);
			drawline.Run(canvasArray[i], propFieldArray[0], omegaArray, sl[i], sr[i]);
			t += dt / (vtkm::Float32)slice + 1.0 / (vtkm::Float32)slice;

		}

		DoNormalize<FieldType,DeviceAdapter> donorm(dim);
		donorm.Run(propFieldArray[0], omegaArray, propFieldArray[1]);


		DoSharpen<FieldType, DeviceAdapter> dosharp(dim);
		dosharp.Run(propFieldArray[1], canvasArray[(loop + 1) % slice]);

		
	}


	saveAs("uflic.pnm", propFieldArray[1], dim[0], dim[1]);
	//vtkm::rendering::Mapper mapper;
	//vtkm::rendering::Canvas canvas(512, 512);
	//vtkm::rendering::Scene scene;

	//scene.AddActor(vtkm::rendering::Actor(
	//	ds.GetCellSet(), ds.GetCoordinateSystem(), ds.GetField(fieldNm), colorTable));
	//vtkm::rendering::Camera camera;
	//SetCamera<ViewType>(camera, ds.GetCoordinateSystem().GetBounds());

	//vtkm::rendering::View2D view;
}
