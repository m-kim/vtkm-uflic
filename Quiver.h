#ifndef QUIVER_H
#define QUIVER_H
#include <vtkm/cont/Algorithm.h>
#include <vtkm/BinaryOperators.h>
#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/VectorAnalysis.h>
#include "Draw.h"

struct Mask
{
  template <typename T, typename U>
  VTKM_EXEC_CONT U operator()(const T& x, const U& y) const
  {
    return vtkm::Dot(x,x) > 0;
  }
};

struct NormalizeFunctor
{
  template <typename T>
  VTKM_EXEC_CONT T operator()(const T& x, const T& y) const
  {
    if (vtkm::Dot(x,x) > 0)
      return vtkm::Normal(x);
    return x;
  }
};

struct Stretch
{
  template <typename T>
  VTKM_EXEC_CONT T operator()(const T&x, const T&y) const
  {
    return x * y;
  }
};

template<typename FieldType, vtkm::Id Size>
class CreateVec : public vtkm::worklet::WorkletMapField
{
public:
  typedef  vtkm::Vec<FieldType, Size> VecType;
  typedef void ControlSignature(FieldIn<>, WholeArrayIn<>, FieldOut<>, FieldOut<>);
  typedef void ExecutionSignature(_1, _2, _3, _4);
  //
  VTKM_CONT
  CreateVec(vtkm::Id &d) : dim(d) {}

  template<typename VecArrayType>
  VTKM_EXEC_CONT void operator()(const VecType pt,
                                 vtkm::Id &idx,
                                 VecArrayType &vecArray,
                                 VecType& out) const
  {
    idx = pt.y * dim + pt.x;
    out = vecArray.Get(idx);
  }

  vtkm::Id dim;
};


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
              CanvasArrayType &canvas,
              vtkm::Id2 spacing)
    {
      using CanvasType = typename std::remove_reference<decltype(canvas)>::type;
      using CanvasValueType = typename CanvasType::ValueType;

      typedef typename vtkm::worklet::DispatcherMapField<DrawLineWorkletType>
        DrawLineWorkletDispatchType;
      vtkm::cont::Algorithm::Transform(
            vec,
            vec,
            vec,
            NormalizeFunctor());

      DrawLineWorkletType drawLineWorklet(dim);
      DrawLineWorkletDispatchType drawLineWorkletDispatch(drawLineWorklet);
      vtkm::cont::ArrayHandleConstant<vtkm::Vec<vtkm::Float32, 2>> length(vtkm::Vec<vtkm::Float32,2>(spacing[0]*0.667, spacing[1] * 0.667), vec.GetNumberOfValues());
      vtkm::cont::Algorithm::Transform(
            vec,
            length,
            vec,
             Stretch());

      typename std::remove_reference<decltype(pl)>::type pr;
      pr.Allocate(pl.GetNumberOfValues());
      vtkm::cont::Algorithm::Transform(pl, vec, pr, vtkm::Sum());

      vtkm::cont::ArrayHandleCounting<vtkm::Int32> idxArray(0,1, pr.GetNumberOfValues());
      vtkm::cont::ArrayHandle<vtkm::Int8> mask;
      vtkm::cont::ArrayHandleConstant<CanvasValueType> origCanvas(255, canvas.GetNumberOfValues());


      mask.Allocate( canvas.GetNumberOfValues());
      vtkm::cont::Algorithm::Transform(
            vec,
            vtkm::cont::ArrayHandleConstant<vtkm::Int8>(0, vec.GetNumberOfValues()),
            mask,  Mask());
      CanvasType omega;
      omega.Allocate(canvas.GetNumberOfValues());
      vtkm::cont::ArrayHandleConstant<CanvasValueType> zero(0, omega.GetNumberOfValues());
      vtkm::cont::ArrayCopy(zero, omega);

      drawLineWorkletDispatch.Invoke(pl, pr, idxArray, mask, canvas, omega, origCanvas);
    }

    DimType dim;
};

#endif
