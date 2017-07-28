#ifndef LucasKanade_H
#define LucasKanade_H

#include <vtkm/Types.h>
#include <vtkm/Math.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/Field.h>
#include <vtkm/cont/ArrayHandle.h>

class LucasKanade : public vtkm::worklet::WorkletMapField
{
public:
  LucasKanade(vtkm::Id2 &_d)
    :dim(_d)
  {

  }

  typedef void ControlSignature(FieldIn<vtkm::Id>,
                                WholeArrayInOut<>,
                                FieldOut<>);


  typedef void ExecutionSignature(_1, _2, _3);


  VTKM_EXEC
  vtkm::Id getIdx(vtkm::Id x, vtkm::Id y) const
  {
    return vtkm::Max(vtkm::Min(y, dim[1]-1),static_cast<vtkm::Id>(0)) * dim[0]
      + vtkm::Max(vtkm::Min(x, dim[0]-1), static_cast<vtkm::Id>(0));
  }

  template<typename WholeArrayInType, typename FieldOutType>
  VTKM_EXEC
  void operator()(const vtkm::Id &idx,
                  const WholeArrayInType &data,
                  FieldOutType &reval
                  ) const
  {
    vtkm::Id x, y;
    y = idx / dim[0];
    x = idx % dim[0];

    FieldOutType up = data[getIdx(x, y + 1)];// / float(ttl);
    FieldOutType right_up = data[getIdx(x + 1, y + 1)];// / float(ttl);
    FieldOutType right = data[getIdx(x + 1, y)];// / float(ttl);
    FieldOutType center = data[getIdx(x, y)];// / float(ttl);


    reval = stencil[0] * center + stencil[1] * up + stencil[2] * right + stencil[3] * right_up;
  }

  vtkm::Id2 dim;
  vtkm::Vec<vtkm::Float32,4> stencil;
};

template<typename FieldType, typename DeviceAdapter>
class DoLucasKanade
{
public:
  DoLucasKanade(vtkm::Id2 _d,
                vtkm::Vec<FieldType, 4> &s)
    : dim(_d),
      stencil(s)
  {

  }

  void Run(
          vtkm::cont::ArrayHandle<FieldType> &in,
          vtkm::cont::ArrayHandle<FieldType> &out) {
    typedef typename vtkm::worklet::DispatcherMapField<LucasKanade>
      LucasKanadeWorkletDispatchType;

    LucasKanade LucasKanadeWorklet(dim, stencil);
    LucasKanadeWorkletDispatchType dispatch(LucasKanadeWorklet);
    in.PrepareForInPlace(DeviceAdapter());
    out.PrepareForInPlace(DeviceAdapter());

    dispatch.Invoke(vtkm::cont::ArrayHandleIndex(dim[0] * dim[1]),
                                    in, out);

  }

  vtkm::Id2 dim;
  vtkm::Vec<FieldType, 4> stencil;
};
#endif
