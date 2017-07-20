//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2014 Sandia Corporation.
//  Copyright 2014 UT-Battelle, LLC.
//  Copyright 2014 Los Alamos National Security.
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================

#ifndef Integrator_h
#define Integrator_h

#include <vtkm/Types.h>


template <typename FieldEvaluateType, typename FieldType, vtkm::IdComponent Size>
class RK4Integrator
{
public:
  VTKM_EXEC_CONT
  RK4Integrator()
    : h(0)
    , h_2(0)
  {
  }

  VTKM_EXEC_CONT
  RK4Integrator(const FieldEvaluateType& field, FieldType _h)
    : f(field)
    , h(_h)
    , h_2(_h / 2.f)
  {
  }

  template <typename PortalType>
  VTKM_EXEC bool Step(const vtkm::Vec<FieldType, Size>& pos,
                      const PortalType& field,
                      vtkm::Vec<FieldType, Size>& out) const
  {
    vtkm::Vec<FieldType, Size> k1, k2, k3, k4;

    if (f.Evaluate(pos, field, k1) && f.Evaluate(pos + h_2 * k1, field, k2) &&
        f.Evaluate(pos + h_2 * k2, field, k3) && f.Evaluate(pos + h * k3, field, k4))
    {
      out = pos + h / 6.0f * (k1 + 2 * k2 + 2 * k3 + k4);
      return true;
    }
    return false;
  }

  FieldEvaluateType f;
  FieldType h, h_2;
};

template <typename FieldEvaluateType, typename FieldType, vtkm::IdComponent Size>
class EulerIntegrator
{
public:
  VTKM_EXEC_CONT
  EulerIntegrator()
    : h(0)
  {
  }
  VTKM_EXEC_CONT
  EulerIntegrator(const FieldEvaluateType& _eval, FieldType _h)
    : eval(_eval)
    , h(_h)
  {
  }

  template <typename PortalType>
  VTKM_EXEC bool Step(const vtkm::Vec<FieldType, Size>& pos,
                      const PortalType& field,
                      vtkm::Vec<FieldType, Size>& out) const
  {
    vtkm::Vec<FieldType, Size> vCur;
//    eval.incrT(h);
    if (eval.Evaluate(pos, field, vCur))
    {
      out = pos + h * vCur;
      return true;
    }
    return false;
  }

  FieldEvaluateType eval;
  FieldType h;
};

#endif // vtk_m_worklet_particleadvection_Integrators_h
