#include <vtkm/Math.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/io/reader/VTKDataSetReader.h>
#include "Evaluator.h"
#include "UFLIC.h"

#include <iostream>
#include <fstream>
#include <chrono>

typedef VTKM_DEFAULT_DEVICE_ADAPTER_TAG DeviceAdapter;

//std::shared_ptr<vtkm::cont::DataSetBuilderExplicitIterative> dataSetBuilder;


int main(int argc, char **argv)
{
  int width, height;

  if (argc < 1){
      width = height = 256;

  }
  else{
      width = height = atoi(argv[1]);
  }


//    treePtr = std::shared_ptr<Tree<DeviceAdapter>>(new Tree<DeviceAdapter>());

//    auto t0 = std::chrono::high_resolution_clock::now();
//    vtkm::rendering::Color bg(0.2f, 0.2f, 0.2f, 1.0f);
//    vtkm::rendering::CanvasRayTracer canvas(width, height);
//    MapperRayTracer mapper(csg.GetCellSet(), cntArray, idxArray, vtxArray );

//    vtkm::rendering::Scene scene;
//    scene.AddActor(vtkm::rendering::Actor(csg.GetCellSet(),
//                                          csg.GetCoordinateSystem(),
//                                          csg.GetField("radius"),
//                                          vtkm::rendering::ColorTable("thermal")));


//    view = new vtkm::rendering::View3D(scene, mapper, canvas, bg);

//    view->Initialize();
//    //glutMainLoop();
//    view->Paint();
//    auto t1 = std::chrono::high_resolution_clock::now();

//    std::cout << "Finished " << width << ": " << std::chrono::duration<double>(t1-t0).count() << "s" << std::endl;
//    view->SaveAs("reg3D.pnm");
}
