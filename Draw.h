#ifndef BRESENHAM_H
#define BRESENHAM_H
#include <vtkm/Types.h>
#include <vtkm/Math.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/Field.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSet.h>
#include "Bounds2.h"

template <typename FieldType, vtkm::IdComponent Size, typename DeviceAdapterTag>
class DrawLine : public vtkm::worklet::WorkletMapField
{
public:
  typedef vtkm::Vec<vtkm::Float32, Size> VecType;

  DrawLine( Bounds &bb,
        vtkm::Id2 d)
    : bounds(bb),
      dim(d)
  {

  }

  typedef void ControlSignature(AtomicArrayInOut<FieldType>,
                                AtomicArrayInOut<FieldType>,
                                FieldIn<VecType>,
                                FieldIn<VecType>,
                                FieldIn<FieldType>);
  typedef void ExecutionSignature(_1, _2, _3, _4, _5);

  template<typename AtomicArrayType>
  void operator()(const AtomicArrayType &canvas,
                  const AtomicArrayType &omega,
                  const VecType& p1,
                  const VecType& p2,
                  FieldType val) const {
    vtkm::Vec<FieldType,Size> p = p1;
    vtkm::Vec<FieldType,Size> d = p2;
    d -= p1;

    float N = vtkm::Max(vtkm::Abs(d.x), vtkm::Abs(d.y));
    if (N < 1e-6){
      N = 1;
    }

    vtkm::Vec<FieldType,Size> s(d.x/N, d.y/N);

    for (int i=0; i<N; i++){
      if (bounds.Contains(vtkm::Round(p))){
        vtkm::Id idx = static_cast<vtkm::Id>(vtkm::Round(p.y))*dim[0] + static_cast<vtkm::Id>(vtkm::Round(p.x));
          canvas.Add(idx, val);//color(255,255,255);
          omega.Add(idx, 1.0);
      }
      p += s;
    }
  }
private:
  Bounds bounds;
  vtkm::Id2 dim;
};



template < typename FieldType, vtkm::IdComponent Size, typename DeviceAdapterTag>
class DrawLineWorklet
{
public:
  typedef vtkm::cont::ArrayHandle<FieldType> FieldHandle;
  typedef FieldType ExecutionTypes<DeviceAdapterTag>::Portal
    FieldPortalType;
  typedef DrawLine< FieldType,
                    Size,
                    DeviceAdapterTag>
    DrawLineWorkletType;

  DrawLineWorklet(vtkm::cont::DataSet &_ds)
    :ds(_ds)
  {

  }


  template <typename PointStorage, typename FieldStorage>
  VTKM_CONT
  void Run(
      vtkm::cont::ArrayHandle<FieldType, FieldStorage> &_canvas,
      vtkm::cont::ArrayHandle<FieldType, FieldStorage> &_omega,

      const vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>, PointStorage>& _pl,
      const vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>, PointStorage>& _pr
    )
  {
    pl = _pl;
    pr = _pr;
    canvas = _canvas.PrepareForInPlace(DeviceAdapterTag());
    omega = _omega.PrepareForInPlace(DeviceAdapterTag());
    run();
  }

  ~DrawLineWorklet() {}

private:
  void run(bool dumpOutput = false)
  {
    typedef typename vtkm::worklet::DispatcherMapField<DrawLineWorkletType>
      DrawLineWorkletDispatchType;

    vtkm::Bounds vounds = ds.GetCoordinateSystem(0).GetBounds();
    Bounds bounds(vounds.X, vounds.Y);
    vtkm::cont::CellSetStructured<2> cells;
    ds.GetCellSet(0).CopyTo(cells);
    vtkm::Id2 dims = cells.GetSchedulingRange(vtkm::TopologyElementTagPoint());


    DrawLineWorkletType drawLineWorklet(bounds, dims);
    DrawLineWorkletDispatchType drawLineWorkletDispatch(drawLineWorklet);
    drawLineWorkletDispatch.Invoke(canvas, omega, pl, pr, canvas);
  }


  vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>> pl, pr;
  vtkm::cont::DataSet ds;
  vtkm::Id maxSteps;
  vtkm::Id ParticlesPerRound;
  FieldPortalType canvas, omega;
};


#endif
