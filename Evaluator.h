#ifndef Evaluator_h
#define Evaluator_h
#include <vtkm/Types.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CellSetStructured.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/DynamicArrayHandle.h>


template <typename PortalType, typename FieldType>
class DoubleGyreField
{
public:
  VTKM_CONT
  DoubleGyreField(const vtkm::Bounds& bb)
    : bounds{ bb },
      omega(2 * PI / 10.0),
      A(0.1),
      epsilon(1e-6),
      t(0.0)
  {
  }

  VTKM_EXEC
  void incrT(FieldType dt){
    t += dt;
  }
  VTKM_EXEC_CONT
  bool Evaluate(const vtkm::Vec<FieldType, 2>& pos,
                const PortalType& vtkmNotUsed(vecData),
                vtkm::Vec<FieldType, 2>& out) const
  {
    if (!bounds.Contains(pos))
      return false;
    out[0] = calcU(pos[0], pos[1], t);
    out[1] = calcV(pos[0], pos[1], t);

    return true;
  }


private:
  VTKM_EXEC
  float a(float t)
  {
   return epsilon * vtkm::Sin(omega * t);
  }
  VTKM_EXEC
  float b(float t)
  {
   return 1 - 2 * epsilon * vtkm::Sin(omega * t);
  }
  VTKM_EXEC
  float f(float x, float t)
  {
    return a(t) * x*x + b(t) * x;
  }
  VTKM_EXEC
  float calcU(float x, float y, float t)
  {
       return -PI * A * vtkm::Sin(PI*f(x,t)) * vtkm::Cos(PI*y);
  }

  VTKM_EXEC
  float calcV(float x, float y, float t)
  {
    return PI * A * vtkm::Cos(PI*f(x,t)) * vtkm::Sin(PI*y) * (2 * a(t) * x + b(t));
  }

  vtkm::Bounds bounds;
  const float omega, A, epsilon;
  float t;

};

#endif
