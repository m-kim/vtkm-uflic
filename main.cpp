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
#include <sstream>

#include "Evaluator.h"
#include "Integrator.h"
#include "UFLIC.h"
#include "Draw.h"
#include "Sharpen.h"
#include "Normalize.h"
#include "Jitter.h"


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
	const vtkm::IdComponent ttl = 4;
  std::vector<vtkm::Vec<VecType, Size>> pl[ttl], pr[ttl];

	for (int i = 0; i < ttl; i++) {
		for (int y = 0; y<dim[0]; y++) {
			for (int x = 0; x<dim[1]; x++) {
				pl[i].push_back(vtkm::Vec<VecType, Size>(x + 0.5, y + 0.5));
				pr[i].push_back(vtkm::Vec<VecType, Size>(x + 0.5, y + 0.5));
			}
		}
	}
	std::vector<vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>>> sl(ttl), sr(ttl);
	for (int i = 0; i<ttl; i++) {
		sl[i] = vtkm::cont::make_ArrayHandle(pl[i]);
		sr[i] = vtkm::cont::make_ArrayHandle(pr[i]);
	}


	vtkm::cont::DataSetBuilderUniform dataSetBuilder;
	vtkm::cont::DataSet ds = dataSetBuilder.Create(dim);

	std::vector<FieldType> canvas[ttl], propertyField[2], omega(dim[0] * dim[1], 0), tex(dim[0] * dim[1], 0);
	vtkm::Float32 t = 0;
	const vtkm::Float32 dt = 0.1;
	
	for (int i = 0; i < 2; i++) {
		propertyField[i].resize(dim[0] * dim[1], 0);
	}

	for (int i = 0; i < ttl; i++) {
		canvas[i].resize(dim[0] * dim[1], 0);
	}
	for (int i = 0; i < canvas[0].size(); i++) {
		tex[i] = canvas[0][i] = rand() % 255;
	}

	vtkm::cont::ArrayHandle<FieldType > canvasArray[ttl], propFieldArray[2], omegaArray, texArray;
	for (int i = 0; i < ttl; i++) {
		canvasArray[i] = vtkm::cont::make_ArrayHandle(&canvas[i][0], canvas[i].size());
	}
	propFieldArray[0] = vtkm::cont::make_ArrayHandle(&propertyField[0][0], propertyField[0].size());
	propFieldArray[1] = vtkm::cont::make_ArrayHandle(&propertyField[1][0], propertyField[1].size());
	omegaArray = vtkm::cont::make_ArrayHandle(&omega[0], omega.size());
	texArray = vtkm::cont::make_ArrayHandle(&tex[0], tex.size());
  vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>> VecArray;
	EvalType eval(t, Bounds(0, dim[0], 0, dim[1]));
	IntegratorType integrator(eval, 3.0);
	ParticleAdvectionWorkletType advect(integrator);
	DrawLineWorkletType drawline(ds);
	DoNormalize<FieldType, DeviceAdapter> donorm(dim);
	DoSharpen<FieldType, DeviceAdapter> dosharp(dim);
	DoJitter<FieldType, DeviceAdapter> dojitter(dim);

	for (int loop = 0; loop < 5; loop++) {
		std::cout << "t: " << t << std::endl;
		for (int i = 0; i < sr[loop % ttl].GetNumberOfValues(); i++) {
			vtkm::Id x, y;
			y = i / dim[0];
			x = i % dim[0];
			sl[loop %ttl].GetPortalControl().Set(i, vtkm::Vec<VecType, 2>(x + 0.5, y + 0.5));
		}
		//reset the current canvas
		for (int i = 0; i < canvasArray[loop % ttl].GetNumberOfValues(); i++) {
			canvasArray[loop % ttl].GetPortalControl().Set(i, rand() % 255);
		}

		for (int i = 0; i < propFieldArray[0].GetNumberOfValues(); i++) {
			propFieldArray[0].GetPortalControl().Set(i, 0);
			propFieldArray[1].GetPortalControl().Set(i, 0);
			omegaArray.GetPortalControl().Set(i, 0);
		}
		for (int i = 0; i < vtkm::Min(ttl, loop+1); i++) {
			advect.Run(sl[i], sr[i], VecArray);
			drawline.Run(canvasArray[i], propFieldArray[0], omegaArray, sl[i], sr[i]);
			//t += dt / (vtkm::Float32)ttl + 1.0 / (vtkm::Float32)ttl;

		}

		sr.swap(sl);

		donorm.Run(propFieldArray[0], omegaArray, propFieldArray[1]);
		std::stringstream fn;
		fn << "uflic-" << loop << ".pnm";
		saveAs(fn.str().c_str(), propFieldArray[1], dim[0], dim[1]);

		//REUSE omegaArray as a temporary cache to sharpen
		dosharp.Run(propFieldArray[1], omegaArray);
		dojitter.Run(omegaArray, texArray, canvasArray[(loop) % ttl]);


	}


	//vtkm::rendering::Mapper mapper;
	//vtkm::rendering::Canvas canvas(512, 512);
	//vtkm::rendering::Scene scene;

	//scene.AddActor(vtkm::rendering::Actor(
	//	ds.GetCellSet(), ds.GetCoordinateSystem(), ds.GetField(fieldNm), colorTable));
	//vtkm::rendering::Camera camera;
	//SetCamera<ViewType>(camera, ds.GetCoordinateSystem().GetBounds());

	//vtkm::rendering::View2D view;
}
