#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/worklet/Normalize.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/rendering/View3D.h>
#include <vtkm/cont/testing/Testing.h>
#include "ReaderUFLIC.h"
#include "MapperUFLIC.h"
#include "ViewUFLIC.h"
#include "CanvasUFLIC.h"
#include "UFLIC.h"
#include "Draw.h"

bool do_print = false;

template< typename VecFld>
vtkm::cont::ArrayHandle<vtkm::Float32> runUFLIC(const vtkm::Id2 &dim,
                                      std::vector<VecFld> &sl,
                                      std::vector<VecFld> &sr,
                                      int ttl = 1,
                                      int loop_cnt = 1
                                      )
{
    const int Size = 2;
    using VecType = vtkm::Float32;
    using FieldType = vtkm::Float32;
    using ArrayType = vtkm::cont::ArrayHandle<FieldType>;


    std::vector<ArrayType> canvasArray;
    ArrayType propFieldArray[2], omegaArray, texArray;

    vtkm::cont::ArrayHandleConstant<vtkm::Id> zero(0, dim[0]*dim[1]);

    propFieldArray[0].Allocate(dim[0] * dim[1]);
    propFieldArray[1].Allocate(dim[0] * dim[1]);
    omegaArray.Allocate(dim[0] * dim[1]);
    texArray.Allocate(dim[0] * dim[1]);

    canvasArray.resize(ttl);
    for (int i = 0; i < ttl; i++) {
      canvasArray[i].Allocate(dim[0] * dim[1]);
    }

    for (int i=0; i<texArray.GetNumberOfValues(); i++){
      texArray.GetPortalControl().Set(i,rand()%255);
    }

    DoNormalize<FieldType> donorm(dim);
    DoSharpen<FieldType> dosharp(dim);
    DoJitter<FieldType> dojitter(dim);
    DrawLineWorklet<FieldType, VecType, Size>  drawline(dim);

    for (int loop = 0; loop < loop_cnt; loop++){
      for (int i=0; i<canvasArray[loop%ttl].GetNumberOfValues(); i++){
        canvasArray[loop%ttl].GetPortalControl().Set(i,rand()%255);
      }
      vtkm::cont::ArrayCopy(zero, propFieldArray[0]);
      vtkm::cont::ArrayCopy(zero, propFieldArray[1]);
      vtkm::cont::ArrayCopy(zero, omegaArray);
      for (int i = 0; i < vtkm::Min(ttl, loop+1); i++) {
        drawline.Run(canvasArray[i], propFieldArray[0], omegaArray, sl[i], sr[i]);
      }

      donorm.Run(propFieldArray[0], omegaArray, propFieldArray[1]);
      if (do_print){
        std::stringstream fn;
        fn << "uflic-" << loop << ".pnm";
        //saveAs(fn.str().c_str(), propFieldArray[1], dim[0], dim[1]);
      }

      //REUSE omegaArray as a temporary cache to sharpen
      dosharp.Run(propFieldArray[1], omegaArray);
      dojitter.Run(omegaArray, texArray, canvasArray[(loop) % ttl]);

    }

    return propFieldArray[0];
}

inline void SetCamera(std::unique_ptr<vtkm::rendering::Camera>& camera,
                      const vtkm::Bounds& coordBounds)
{
  camera = std::unique_ptr<vtkm::rendering::Camera>(new vtkm::rendering::Camera());
  camera->ResetToBounds(coordBounds);
  camera->Azimuth(static_cast<vtkm::Float32>(45.0));
  camera->Elevation(static_cast<vtkm::Float32>(45.0));
}


void Render(ViewUFLIC& view)
{
  view.Initialize();
  view.Paint();
}
inline vtkm::cont::DataSet readVTKDataSet(const char* fname)
{
  vtkm::cont::DataSet ds;
  vtkm::io::reader::VTKDataSetReader reader(fname);
  try
  {
    ds = reader.ReadDataSet();
  }
  catch (vtkm::io::ErrorIO& e)
  {
    std::string message("Error reading: ");
    message += fname;
    message += ", ";
    message += e.GetMessage();

    VTKM_TEST_FAIL(message.c_str());
  }

  return ds;
}
void addField(vtkm::cont::DataSet &ds){
    auto cnt = ds.GetField(0).GetData().GetNumberOfValues();
    vtkm::cont::DataSetFieldAdd dsf;
    std::vector<vtkm::Float32> vars(cnt);
    for (int i=0; i<cnt ;i++){
        vars[i] = vtkm::Float32(i)/vtkm::Float32(cnt);
    }


    dsf.AddPointField(ds, "pointvar", vars);

}

vtkm::cont::ColorTable buildColorTableFromArray(const vtkm::cont::DynamicArrayHandle &_h){

    auto handle = _h.Cast<vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,3>>>();
    vtkm::worklet::DispatcherMapField<vtkm::worklet::Normalize>(
        vtkm::worklet::Normalize())
        .Invoke(handle);

    std::vector<double> out(_h.GetNumberOfValues()*3);
    for (int i=0; i<handle.GetNumberOfValues(); i++){
        auto tmp = handle.GetPortalConstControl().Get(i);
        out[i*3] = tmp[0];
        out[i*3+1] = tmp[1];
        out[i*3+2] = tmp[2];
    }


    std::string name = "inplace";
    vtkm::Vec<double,3> nanIs(0,0,0);
     return vtkm::cont::ColorTable(name, vtkm::cont::ColorSpace::RGB, nanIs,
                              out);

}

auto mapper = std::unique_ptr<MapperUFLIC>();
auto canvas = std::unique_ptr<CanvasUFLIC>();
auto scene = std::unique_ptr<vtkm::rendering::Scene>();
auto camera = std::unique_ptr<vtkm::rendering::Camera>();
auto view = std::unique_ptr<ViewUFLIC>();

int main(){
  const vtkm::Id2 dim(512,512);
  mapper = std::unique_ptr<MapperUFLIC>(new MapperUFLIC());
  canvas = std::unique_ptr<CanvasUFLIC>(new CanvasUFLIC(dim[0],dim[1]));
  scene = std::unique_ptr<vtkm::rendering::Scene>(new vtkm::rendering::Scene());
  camera = std::unique_ptr<vtkm::rendering::Camera>(new vtkm::rendering::Camera());
  view = std::unique_ptr<ViewUFLIC>(new ViewUFLIC(*scene, *mapper, *canvas, *camera, vtkm::rendering::Color(0,0,0,1), vtkm::rendering::Color(1,0,0,1)));

  mapper->SetShadingOn(false);
  auto ds = readVTKDataSet("edelta-velocity.vtk");
  addField(ds);
  static std::string fieldNm = "pointvar";

  scene->AddActor(vtkm::rendering::Actor(
  ds.GetCellSet(), ds.GetCoordinateSystem(), ds.GetField(fieldNm), buildColorTableFromArray(ds.GetField("vectors").GetData())));
  SetCamera(camera, ds.GetCoordinateSystem().GetBounds());
  view = std::unique_ptr<ViewUFLIC>(new ViewUFLIC(*scene, *mapper, *canvas, *camera, vtkm::rendering::Color(0,0,0,1), vtkm::rendering::Color(1,0,0,1)));
  Render(*view);

  std::vector<vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 2>>> sl, sr;
  std::vector<vtkm::cont::ArrayHandle<vtkm::Id>> pIdx;
  sr.push_back(canvas->pixelPos);
  sl.push_back(canvas->pixelPrePos);
  sr[0].GetPortalConstControl().Get(0);
  sl[0].GetPortalConstControl().Get(0);
  pIdx.push_back(canvas->pixelIdx);

  runUFLIC(dim, sl, sr);

return 0;
}
