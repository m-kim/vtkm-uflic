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
#include "ScreenSpaceLIC.h"
#include "Draw.h"

bool do_print = true;

template< typename VecFld>
vtkm::cont::ArrayHandle<vtkm::Int32> draw(const vtkm::Id2 &dim,
                                      std::vector<VecFld> &sl,
                                      std::vector<VecFld> &sr,
                                      vtkm::cont::ArrayHandle<vtkm::Float32> &depth,
                                          vtkm::Float32 stepsize = 1.0,
                                      int ttl = 1,
                                      int loop_cnt = 1
                                      )
{
    ScreenSpaceLIC<VectorField< vtkm::Float32,2>, vtkm::Float32> lic(dim, stepsize, ttl, loop_cnt);
    return lic.draw(sl, sr, depth);
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

  sr.push_back(canvas->pixelPos);
  sl.push_back(canvas->pixelPrePos);

  draw(dim, sl, sr, canvas->GetDepthBuffer(), 0.03);

return 0;
}
