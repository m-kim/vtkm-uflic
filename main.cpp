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
#include <memory>

#include "UFLIC.h"

#include "Reader.h"

#include <iostream>
#include <fstream>
#include <chrono>
#include <type_traits>

std::tuple<int,int,int, bool>
parse(int argc, char **argv){
  bool do_print = false;
  int x = 512;
  int y = 256;
  int which = 0;

  for (int i=1; i<argc; i++){
    if (!strcmp(argv[i], "bfield")){
      which = 1;
    }
    else if (!strcmp(argv[i], "PSI")){
      which = 2;
    }
    else if (!strcmp(argv[i], "dims")){
      if (i+1 < argc && i+2 < argc){
        x = atoi(argv[i+1]);
        y = atoi(argv[i+2]);
        i += 2;
      }
    }
    else if (!strcmp(argv[i], "print")){
      do_print = true;
    }
  }

  return std::make_tuple(which,x,y, do_print);
}

int main(int argc, char **argv)
{
  const int Size = 2;
  typedef vtkm::Float32 VecType;

  auto ret = parse(argc,argv);

  std::shared_ptr<Reader<VecType, Size>> reader;

  if (std::get<0>(ret) == 1){
    reader = std::shared_ptr<ReaderVTK<VecType, Size>>(new ReaderVTK<VecType, Size>("/home/ybk/Projects/vtkm-uflic/BField_2d.vtk", 12));
    UFLIC<VectorField<VecType,Size>,VecType,Size> uflic;
    uflic.do_print = std::get<3>(ret);
    uflic.run(reader);
  }

  else if (std::get<0>(ret) == 2){
    //std::shared_ptr<Reader<VecType, Size,  ReaderPS<VecType, Size,ReaderXGC<VecType,Size>>>> reader(new ReaderPS<VecType, Size, ReaderXGC<VecType,Size>>("/home/mkim/vtkm-uflic/psi2q/2D_packed/psi2D_packed_normalized_256_99.vec", vtkm::Id2(256,256), Bounds(0,256,0,256)));
    reader = std::shared_ptr<ReaderXGC<VecType, Size>>(new ReaderXGC<VecType, Size>("/home/ybk/Projects/vtkm-uflic/psi2q/2D_packed/psi2D_packed_512_", vtkm::Id2(512,512), Bounds(0,512,0,512), 12));
    UFLIC<VectorField<VecType,Size>,VecType,Size> uflic;
    uflic.do_print = std::get<3>(ret);
    uflic.run(reader);
  }
  else{

    int x = std::get<1>(ret);
    int y = std::get<2>(ret);
    reader = std::shared_ptr<ReaderCalc<VecType, Size>>(new ReaderCalc<VecType, Size>("XGC_", vtkm::Id2(x,y), Bounds(0,x,0,y), vtkm::Vec<VecType,Size>(2,1), 12));
    UFLIC<DoubleGyreField<VecType,Size>,VecType,Size> uflic;
    uflic.do_print = std::get<3>(ret);
    uflic.run(reader);
  }
}
