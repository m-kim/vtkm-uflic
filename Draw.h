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
                                FieldIn<>,
                                WholeArrayIn<>,
#if 0
                                AtomicArrayInOut<>,
                                AtomicArrayInOut<>,
#else
                                WholeArrayInOut<>,
                                WholeArrayInOut<>,
#endif

                                WholeArrayIn<>);
  typedef void ExecutionSignature(_1, _2, _3, _4, _5, _6, _7);

  template<typename AtomicArrayType,
           typename MaskArrayType,
           typename IdxType,
           typename PrevCanvasArray>
  VTKM_EXEC
  void operator()(
                  const VecType& p1,
                  const VecType& p2,
                  const IdxType &idxArray,
                  MaskArrayType &mask,
                  AtomicArrayType &canvas,
                  AtomicArrayType &omega,
                  const PrevCanvasArray &valarray) const {
    if (!outside(p1) && !outside(p2)){

			vtkm::Vec<VecComponentType, Size> p = p1;
      const vtkm::Vec<VecComponentType, Size> d = (p2 - p);
      if (vtkm::Magnitude(d) > 1e-6){
        float N = vtkm::Max(vtkm::Abs(d[0]), vtkm::Abs(d[1]));
        if (N > 0) {
          auto val = valarray[idxArray];
          const vtkm::Vec<VecComponentType, Size> s(d[0] / N, d[1] / N);

          vtkm::Id idx = static_cast<vtkm::Id>(vtkm::Round(p[1]))*dim[0] + static_cast<vtkm::Id>(vtkm::Round(p[0]));
          for (int i = 0; i<N && mask[idx] > 0; i++) {
            if (!outside(vtkm::Round(p))) {
    #if 0
              canvas.Add(idx, val);//color(255,255,255);
              omega.Add(idx, 1);
    #else
              auto cnvs = canvas.Get(idx);
              cnvs += val;
              canvas.Set(idx, cnvs);//color(255,255,255);
              auto newidx = omega.Get(idx);

              omega.Set(idx, newidx + 1);
    #endif
            }
            p += s;
            idx = static_cast<vtkm::Id>(vtkm::Round(p[1]))*dim[0] + static_cast<vtkm::Id>(vtkm::Round(p[0]));
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


  template <typename PointStorage,
            typename FieldStorage,
            typename MaskArrayHandleType,
            typename IdxArrayType>
  VTKM_CONT
  void Run(
      vtkm::cont::ArrayHandle<FieldType, FieldStorage> &_canvas0,
			vtkm::cont::ArrayHandle<FieldType, FieldStorage> &_canvas1,
      vtkm::cont::ArrayHandle<FieldType, FieldStorage> &_omega,
      MaskArrayHandleType &_mask,
      IdxArrayType &_idxArray,
      const vtkm::cont::ArrayHandle<vtkm::Vec<VecComponentType, Size>, PointStorage>& _pl,
      const vtkm::cont::ArrayHandle<vtkm::Vec<VecComponentType, Size>, PointStorage>& _pr
    )
  {


    typedef typename vtkm::worklet::DispatcherMapField<DrawLineWorkletType>
      DrawLineWorkletDispatchType;

    DrawLineWorkletType drawLineWorklet(dims);
    DrawLineWorkletDispatchType drawLineWorkletDispatch(drawLineWorklet);
    drawLineWorkletDispatch.Invoke(_pl, _pr, _idxArray, _mask, _canvas1, _omega, _canvas0);
  }

  ~DrawLineWorklet() {}

private:



  vtkm::Id2 dims;
  vtkm::Id maxSteps;
  vtkm::Id ParticlesPerRound;
};


#endif
