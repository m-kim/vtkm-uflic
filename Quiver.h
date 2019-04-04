#ifndef QUIVER_H
#define QUIVER_H
#include <vtkm/cont/Algorithm.h>
#include <vtkm/BinaryOperators.h>
#include "Draw.h"

template<
         typename VecArrayType
         ,typename CanvasArrayType
         ,typename DimType>
class Quiver
{
public:
  using CanvasType = typename CanvasArrayType::ValueType;
  using VecType = typename VecArrayType::ValueType::ComponentType;


  typedef DrawLine< CanvasType,
                    VecType,
                    2>
    DrawLineWorkletType;

    Quiver(DimType d)
      :dim(d)
    {}
    ~Quiver(){}

    void draw(VecArrayType &pl,
              VecArrayType &vec,
              CanvasArrayType &canvas)
    {

      typedef typename vtkm::worklet::DispatcherMapField<DrawLineWorkletType>
        DrawLineWorkletDispatchType;

      DrawLineWorkletType drawLineWorklet(dim);
      DrawLineWorkletDispatchType drawLineWorkletDispatch(drawLineWorklet);

      typename std::remove_reference<decltype(pl)>::type pr;
      pr.Allocate(pl.GetNumberOfValues());
      vtkm::cont::Algorithm::Transform(pl, vec, pl, vtkm::Sum());

      vtkm::cont::ArrayHandleCounting<vtkm::Int32> idxArray(0,1, pr.GetNumberOfValues());
      vtkm::cont::ArrayHandleConstant<vtkm::Int8> mask(1, pr.GetNumberOfValues());
      vtkm::cont::ArrayHandleConstant<CanvasType> origCanvas(255, pr.GetNumberOfValues());

      typename std::remove_reference<decltype(canvas)>::type omega;
      omega.Allocate(pl.GetNumberOfValues());
      drawLineWorkletDispatch.Invoke(pl, pr, idxArray, mask, canvas, omega, origCanvas);
    }

    DimType dim;
};

#endif
