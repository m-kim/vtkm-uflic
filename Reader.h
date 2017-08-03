#ifndef READER_H
#define READER_H
#include <iostream>
#include <sstream>
#include "Bounds2.h"
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetBuilderUniform.h>

template <typename VecType, vtkm::Id Size, class Derived>
class Reader
{
public:
  typedef Derived DerivedType;
  typedef VectorField<VecType> EvalType;

  Reader(std::string fn,
         vtkm::Id2 d,
         Bounds bb,
         vtkm::Vec<VecType, Size> sp)
    : filename(fn),
      dim(d),
      bounds(bb),
        spacing(sp)
  {

  }
  virtual void read(std::vector<vtkm::Vec<VecType, Size>> &in) = 0;

  virtual void next(std::vector<vtkm::Vec<VecType, Size>> &in){}
  vtkm::Id2 dim;
  Bounds bounds;
  const vtkm::Vec<VecType, Size> spacing;
  vtkm::cont::DataSetBuilderUniform dataSetBuilder;
  vtkm::cont::DataSet ds;
  std::stringstream filename;
};

template <typename VecType, vtkm::Id Size, class Derived>
class ReaderPS : public Reader<VecType, Size, ReaderPS<VecType, Size, Derived>>
{
public:
  ReaderPS(std::string fn,
           vtkm::Id2 d,
           Bounds bb)
    : Reader<VecType, Size, ReaderPS>(fn,
                            d,
                            bb,
                            vtkm::Vec<VecType, Size>(1,1)
                            )
  {

    this->ds = this->dataSetBuilder.Create(this->dim);
  }

  void read(std::vector<vtkm::Vec<VecType, Size>> &in)
  {
    std::cout << this->filename.str() << std::endl;

    std::string line;
    std::ifstream file(this->filename.str());

    while (std::getline(file, line)) {
      std::stringstream ss;
      ss << line;
      std::string tok;
      vtkm::Vec<VecType, Size> vec;
      //while (std::getline(ss, tok, ' ')) {
      std::getline(ss, tok, ' ');
      vec[0] = atof(tok.c_str());
      std::getline(ss, tok, ' ');
      vec[1] = atof(tok.c_str());

      in.push_back(vec);
      //}
    }
    //	String text = null;

    //		String[] subtext = splitTokens(text, " ");

    //		vecs[cnt].x = float(subtext[0]);
    //		vecs[cnt].y = float(subtext[1]);
    //		cnt += 1;
    //	}
    //}
    //catch (IOException e) {
    //	e.printStackTrace();
    //}
  }
};

template <typename VecType, vtkm::Id Size>
class ReaderVTK : public Reader<VecType, Size, ReaderVTK<VecType, Size>>
{
public:
  ReaderVTK(std::string fn)
    : Reader<VecType, Size, ReaderVTK>(fn,
                            vtkm::Id2(512,512),
                            Bounds(0,512, 0, 512),
                            vtkm::Vec<VecType, Size>(1,1)
                            )
  {

  }
  void read(std::vector<vtkm::Vec<VecType, Size>> &in)
  {
    std::cout << this->filename.str() << std::endl;
    ds = readVTK(this->filename.str());
  //  vtkm::Bounds vounds = ds.GetCoordinateSystem(0).GetBounds();
  //  Bounds bounds(vounds.X, vounds.Y);

    vtkm::cont::CellSetStructured<3> cells;
    ds.GetCellSet(0).CopyTo(cells);
    dim3 = cells.GetSchedulingRange(vtkm::TopologyElementTagPoint());
    this->dim = vtkm::Id2(dim3[0], dim3[2]);
    vtkm::cont::Field fld = ds.GetField(0);
    dah = fld.GetData();
    ah = dah.Cast<vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64,3>>>();
    std::cout << ah.GetNumberOfValues() << std::endl;
    loop = 0;
    next(in);
  }

  void next(std::vector<vtkm::Vec<VecType, Size>> &in)
  {
    in.resize(0);
    for (int z=0; z<this->dim[1]; z++){
      for (int y=loop; y<loop+1; y++){
        for (int x=0; x<this->dim[0]; x++){
          vtkm::Id idx = z*this->dim[0]*128 + y*this->dim[0]+x;
          vtkm::Vec<vtkm::Float64,3> vec  = ah.GetPortalConstControl().Get(idx);
          in.push_back(vtkm::Vec<vtkm::Float32, 2>(vec[0],vec[2]));
        }
      }
    }
    loop++;

  }

  vtkm::cont::DataSet readVTK(std::string fn)
  {
    vtkm::cont::DataSet ds;
    vtkm::io::reader::VTKDataSetReader rdr(fn.c_str());
    try
    {
      ds = rdr.ReadDataSet();
    }
    catch (vtkm::io::ErrorIO &e) {
      std::string message("Error reading: ");
      message += fn.c_str();
      message += ", ";
      message += e.GetMessage();
      std::cerr << message << std::endl;
    }
    return ds;

  }


  vtkm::cont::DataSet ds;
  vtkm::cont::DynamicArrayHandle dah;
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64,3>> ah;

  vtkm::Id3 dim3;
  vtkm::Id loop;
};

template <typename VecType, vtkm::Id Size>
class ReaderCalc : public Reader<VecType, Size, ReaderCalc<VecType, Size>>
{
public:
  typedef DoubleGyreField<VecType, Size> EvalType;
  ReaderCalc(std::string fn,
             vtkm::Id2 d = vtkm::Id2(512,512),
             Bounds bb = Bounds(0,512, 0, 512),
             vtkm::Vec<VecType, Size> sp = vtkm::Vec<VecType, Size>(2,1))
    : Reader<VecType, Size, ReaderCalc>(fn,
                            d,
                            bb,
                            sp
                            )
  {
  }
  void read(std::vector<vtkm::Vec<VecType, Size>> &in)
  {
  }
};


template <typename VecType, vtkm::Id Size>
class ReaderXGC : public ReaderPS<VecType, Size, ReaderXGC<VecType, Size>>
{
public:
  ReaderXGC(std::string fn,
           vtkm::Id2 d,
           Bounds bb)
    : ReaderPS<VecType, Size,ReaderXGC>(fn,
                            d,
                            bb
                            ),
      base_fn(fn)

  {
    this->ds = this->dataSetBuilder.Create(this->dim);
    loop = 3;
    this->filename << base_fn << loop << ".vel";
    std::cout << this->filename.str() << std::endl;
    loop++;
  }

  void next(std::vector<vtkm::Vec<VecType, Size>> &in)
  {
    in.resize(0);
    this->filename.str("");
    this->filename << base_fn << loop << ".vel";
    this->read(in);
    loop++;
  }

  std::string base_fn;
  int loop;
};
#endif
