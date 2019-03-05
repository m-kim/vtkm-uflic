#ifndef READERUFLIC_H
#define READERUFLIC_H
#include <vtkm/Types.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/rendering/Canvas.h>
#include "Evaluator.h"
#include "Reader.h"

template <typename VecType, vtkm::Id Size>
class ReaderUFLIC : public Reader<VecType, Size>
{
public:
  using ColorBufferType = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 4>>;
  using EvalType = VectorField<VecType>;

  using ArrayType = vtkm::cont::ArrayHandle<vtkm::Vec<VecType, Size>>;
  
  ReaderUFLIC(vtkm::rendering::Canvas& _canvas)
    : Reader<VecType, Size>("",
                            vtkm::Id2(_canvas.GetWidth(),_canvas.GetHeight()),
                            Bounds(0, _canvas.GetWidth(), 0, _canvas.GetWidth()),
                            vtkm::Vec<VecType, Size>(1,1),
                            10
                            )
    ,canvas(_canvas)
  {
      std::cout << this->dim << std::endl;
  }
  void readFile()
  {
  
    this->ah.Allocate(canvas.GetColorBuffer().GetNumberOfValues());
    std::cout << this->ah.GetNumberOfValues() << std::endl;
    parse(this->ah);

  }

  void next(ArrayType  &in)
  {
    in = this->ah;
  }

  void parse(ArrayType in)
  {
    ColorBufferType::PortalConstControl colorPortal = canvas.GetColorBuffer().GetPortalConstControl();
    std::cout << this->dim << std::endl;  
    for (int y=0; y<this->dim[1]; y++){
        for (int x=0; x<this->dim[0]; x++){
            int idx = y*this->dim[0] + x;
              auto vec  = colorPortal.Get(idx);
              in.GetPortalControl().Set(idx, vtkm::Vec<vtkm::Float32, 2>(vec[0],vec[1]));
        }
        
    }
    }
    vtkm::rendering::Canvas canvas;
    ArrayType ah;
};

#endif
