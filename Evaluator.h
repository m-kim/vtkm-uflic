#ifndef Evaluator_h
#define Evaluator_h
#include <vtkm/Types.h>
#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/CellSetStructured.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/DynamicArrayHandle.h>
#include "Bounds2.h"


template <typename PortalType, typename FieldType>
class DoubleGyreField
{
public:
  VTKM_CONT
  DoubleGyreField()
    : omega(2 * vtkm::Pi() / 10.0),
      A(0.1),
      epsilon(1e-6),
      t(0.0)
  {
  }

  VTKM_CONT
  DoubleGyreField(const Bounds& bb)
    : bounds( bb ),
      omega(2 * vtkm::Pi() / 10.0),
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
    out[0] = pos[0] + calcU(pos[0], pos[1], t);
    out[1] = pos[1] + calcV(pos[0], pos[1], t);

    return true;
  }


private:
  VTKM_EXEC
  FieldType a(FieldType t) const
  {
   return epsilon * vtkm::Sin(omega * t);
  }
  VTKM_EXEC
  FieldType b(FieldType t) const
  {
   return 1 - 2 * epsilon * vtkm::Sin(omega * t);
  }
  VTKM_EXEC
  FieldType f(FieldType x, FieldType t) const
  {
    return a(t) * x*x + b(t) * x;
  }
  VTKM_EXEC
  FieldType calcU(FieldType x, FieldType y, FieldType t) const
  {
       return -vtkm::Pi() * A * vtkm::Sin(vtkm::Pi()*f(x,t)) * vtkm::Cos(vtkm::Pi()*y);
  }

  VTKM_EXEC
  FieldType calcV(FieldType x, FieldType y, FieldType t) const
  {
    return vtkm::Pi() * A * vtkm::Cos(vtkm::Pi()*f(x,t)) * vtkm::Sin(vtkm::Pi()*y) * (2 * a(t) * x + b(t));
  }

  Bounds bounds;
  FieldType omega, A, epsilon;
  FieldType t;

};

#endif
