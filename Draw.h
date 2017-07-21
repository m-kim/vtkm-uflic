#ifndef BRESENHAM_H
#define BRESENHAM_H
#include <vtkm/Types.h>
#include <vtkm/Math.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/Field.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSet.h>
#include "Bounds2.h"

template <typename FieldType, typename VecComponentType, vtkm::IdComponent Size, typename DeviceAdapterTag>
class DrawLine : public vtkm::worklet::WorkletMapField
{
public:
  typedef vtkm::Vec<VecComponentType, Size> VecType;

  DrawLine( Bounds &bb,
        vtkm::Id2 d)
    : bounds(bb),
      dim(d)
  {

  }

#if 1
  typedef void ControlSignature(AtomicArrayInOut<FieldType>,
                                AtomicArrayInOut<FieldType>,
#else
	typedef void ControlSignature(WholeArrayInOut<FieldType>,
																WholeArrayInOut<FieldType>,
#endif
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
    vtkm::Vec<VecComponentType,Size> p = p1;
    const vtkm::Vec<VecComponentType,Size> d = p2 - p;

    float N = vtkm::Max(vtkm::Abs(d[0]), vtkm::Abs(d[1]));
    if (N < 1e-6){
      N = 1;
    }

    const vtkm::Vec<VecComponentType,Size> s(d[0]/N, d[1]/N);

    for (int i=0; i<N; i++){
      if (bounds.Contains(vtkm::Round(p))){
        vtkm::Id idx = static_cast<vtkm::Id>(vtkm::Round(p[1]))*dim[0] + static_cast<vtkm::Id>(vtkm::Round(p[0]));
#if 1
				canvas.Add(idx, val);//color(255,255,255);
				omega.Add(idx, 1);
#else
				canvas.Set(idx, val);//color(255,255,255);
				omega.Set(idx, 1);
#endif
			}
      p += s;
    }
	}
private:
  Bounds bounds;
  vtkm::Id2 dim;
};



template < typename FieldType, typename VecComponentType, vtkm::IdComponent Size, typename DeviceAdapterTag>
class DrawLineWorklet
{
public:
	typedef DrawLine< FieldType,
										VecComponentType,
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
      vtkm::cont::ArrayHandle<FieldType, FieldStorage> &_canvas0,
			vtkm::cont::ArrayHandle<FieldType, FieldStorage> &_canvas1,
      vtkm::cont::ArrayHandle<FieldType, FieldStorage> &_omega,

      const vtkm::cont::ArrayHandle<vtkm::Vec<VecComponentType, Size>, PointStorage>& _pl,
      const vtkm::cont::ArrayHandle<vtkm::Vec<VecComponentType, Size>, PointStorage>& _pr
    )
  {
    pl = _pl;
    pr = _pr;
		canvas[0] = _canvas0;
		canvas[1] = _canvas1;
		omega = _omega;

		canvas[0].PrepareForInPlace(DeviceAdapterTag());
		canvas[1].PrepareForInPlace(DeviceAdapterTag());
		omega.PrepareForInPlace(DeviceAdapterTag());
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
    drawLineWorkletDispatch.Invoke(canvas[1], omega, pl, pr, canvas[0]);
  }


  vtkm::cont::ArrayHandle<vtkm::Vec<VecComponentType, Size>> pl, pr;
  vtkm::cont::DataSet ds;
  vtkm::Id maxSteps;
  vtkm::Id ParticlesPerRound;
  vtkm::cont::ArrayHandle<FieldType> canvas[2], omega;
};


#endif
