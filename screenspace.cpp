#include <vtkm/rendering/Actor.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/worklet/Normalize.h>
#include <vtkm/cont/DataSetBuilderExplicit.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include <vtkm/rendering/View3D.h>
#include <vtkm/cont/testing/Testing.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include "ReaderUFLIC.h"
#include "MapperUFLIC.h"
#include "ViewUFLIC.h"
#include "CanvasUFLIC.h"
#include "ScreenSpaceLIC.h"
#include "Draw.h"

bool do_print = true;

template< typename VecFld, typename VecIdx>
void draw(const vtkm::Id2 &dim,
          VecFld &vecArray,
          VecIdx &indexArray,
          vtkm::cont::ArrayHandle<vtkm::Float32> &depth,
              vtkm::Float32 stepsize = 1.0,
          int ttl = 1,
          int loop_cnt = 1
          )
{
    ScreenSpaceLIC<VectorField< vtkm::Float32,2>, vtkm::Float32> lic(dim, stepsize, ttl, loop_cnt);
    lic.draw(vecArray, indexArray, depth);
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
  auto ds = readVTKDataSet("/home/mark/external/Dropbox/Presentation/ice-train-vel.vtk");
  addField(ds);
  static std::string fieldNm = "pointvar";

  scene->AddActor(vtkm::rendering::Actor(
  ds.GetCellSet(), ds.GetCoordinateSystem(), ds.GetField(fieldNm), buildColorTableFromArray(ds.GetField("vectors").GetData())));
  SetCamera(camera, ds.GetCoordinateSystem().GetBounds());
  view = std::unique_ptr<ViewUFLIC>(new ViewUFLIC(*scene, *mapper, *canvas, *camera, vtkm::rendering::Color(0,0,0,1), vtkm::rendering::Color(1,0,0,1)));
  Render(*view);

  vtkm::cont::ArrayHandle<vtkm::Id> conn;
  conn.Allocate(dim[0]*dim[1]);
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 3>> velField;
  velField.Allocate(dim[0]*dim[1]);
  for (int i=0; i<canvas->pixelVel.GetNumberOfValues(); i++){
    auto val2 = canvas->pixelVel.GetPortalConstControl().Get(i);
   vtkm::Vec<vtkm::Float32, 3> fin2(val2[0], val2[1], 0.0);;
    velField.GetPortalControl().Set(i, fin2);
    conn.GetPortalControl().Set(i, i);
  }

  vtkm::cont::DataSetBuilderExplicit dsb;
  vtkm::cont::DataSet outds = dsb.Create(velField, vtkm::CellShapeTagVertex(), 1, conn);

  vtkm::io::writer::VTKDataSetWriter writer("output.vtk");
  writer.WriteDataSet(outds);
//  draw(dim, canvas->pixelVel, canvas->pixelIdx, canvas->GetDepthBuffer(), 4, 10, 10);

return 0;
}
