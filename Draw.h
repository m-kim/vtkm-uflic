#ifndef BRESENHAM_H
#define BRESENHAM_H
#include <vtkm/Types.h>
#include <vtkm/Math.h>
#include <vtkm/worklet/WorkletMapField.h>
#include <vtkm/cont/Field.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/DataSet.h>

template <typename FieldType, typename VecComponentType, vtkm::IdComponent Size>
class DrawLine : public vtkm::worklet::WorkletMapField
{
public:
  typedef vtkm::Vec<VecComponentType, Size> VecType;

  DrawLine(
        vtkm::Id2 d)
    :
      dim(d)
  {

  }

  VTKM_EXEC
  bool outside(const VecType &pos) const
  {
    return pos[0] < 0 || pos[0] >= dim[0] || pos[1] < 0 || pos[1] >= dim[1] || pos[0] != pos[0] || pos[1] != pos[1];
  }
  typedef void ControlSignature(
                                FieldIn<>,
                                FieldIn<>,
#if 0
                                AtomicArrayInOut<>,
                                AtomicArrayInOut<>,
#else
                                WholeArrayInOut<>,
                                WholeArrayInOut<>,
#endif
                                WholeArrayInOut<>);
  typedef void ExecutionSignature(_1, _2, _3, _4, _5, WorkIndex);

  template<typename AtomicArrayType>
  VTKM_EXEC
  void operator()(
                  const VecType& p1,
                  const VecType& p2,
                  AtomicArrayType &canvas,
                  AtomicArrayType &omega,
                  const AtomicArrayType &valarray,
                  const vtkm::Id &widx) const {
    if (!outside(p1) && !outside(p2)){

			vtkm::Vec<VecComponentType, Size> p = p1;
			const vtkm::Vec<VecComponentType, Size> d = p2 - p;
      if (vtkm::Magnitude(d) > 1e-6){
        float N = vtkm::Max(vtkm::Abs(d[0]), vtkm::Abs(d[1]));
        if (N > 0) {
          auto val = valarray[widx];
          const vtkm::Vec<VecComponentType, Size> s(d[0] / N, d[1] / N);

          for (int i = 0; i<N; i++) {
            if (!outside(vtkm::Round(p))) {
              const vtkm::Id idx = static_cast<vtkm::Id>(vtkm::Round(p[1]))*dim[0] + static_cast<vtkm::Id>(vtkm::Round(p[0]));
    #if 0
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

      }
    }
	}
private:
  const vtkm::Id2 dim;
};



template < typename FieldType, typename VecComponentType, vtkm::IdComponent Size>
class DrawLineWorklet
{
public:
	typedef DrawLine< FieldType,
										VecComponentType,
                    Size>
    DrawLineWorkletType;

  DrawLineWorklet(
                  vtkm::Id2 d)
    :
      dims(d)
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

    //canvas[0].PrepareForInPlace(DeviceAdapterTag());
    //canvas[1].PrepareForInPlace(DeviceAdapterTag());
    //omega.PrepareForInPlace(DeviceAdapterTag());
    run();
  }

  ~DrawLineWorklet() {}

private:
  void run(bool dumpOutput = false)
  {
    typedef typename vtkm::worklet::DispatcherMapField<DrawLineWorkletType>
      DrawLineWorkletDispatchType;

    DrawLineWorkletType drawLineWorklet(dims);
    DrawLineWorkletDispatchType drawLineWorkletDispatch(drawLineWorklet);
    drawLineWorkletDispatch.Invoke(pl, pr,canvas[1], omega, canvas[0]);
  }


  vtkm::cont::ArrayHandle<vtkm::Vec<VecComponentType, Size>> pl, pr;

  vtkm::Id2 dims;
  vtkm::Id maxSteps;
  vtkm::Id ParticlesPerRound;
  vtkm::cont::ArrayHandle<FieldType> canvas[2], omega;
};


#endif
