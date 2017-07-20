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

#ifndef vtk_m_worklet_particleadvection_ParticleAdvectionWorklets_h
#define vtk_m_worklet_particleadvection_ParticleAdvectionWorklets_h

#include <vtkm/Types.h>
#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleCounting.h>
#include <vtkm/cont/CellSetExplicit.h>
#include <vtkm/cont/CellSetStructured.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/Field.h>
#include <vtkm/exec/ExecutionObjectBase.h>

#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/particleadvection/Particles.h>



template <typename IntegratorType, typename FieldType, vtkm::IdComponent Size, typename DeviceAdapterTag>
class ParticleAdvectWorklet : public vtkm::worklet::WorkletMapField
{
public:
  typedef vtkm::Vec<vtkm::Float32, Size> VecType;
  typedef vtkm::cont::ArrayHandle<VecType> FieldHandle;
  typedef typename FieldHandle::template ExecutionTypes<DeviceAdapterTag>::PortalConst
    FieldPortalConstType;

  typedef void ControlSignature(FieldIn<VecType>, FieldOut<VecType>);
  typedef void ExecutionSignature(_1, _2);
  typedef _1 InputDomain;

  VTKM_EXEC void operator()(const VecType& p1, VecType& p2) const
  {

    integrator.Step(p1, field, p2);
  }

  ParticleAdvectWorklet(const IntegratorType& it, const FieldPortalConstType& f)
    : integrator(it)
    , field(f)
  {
  }

  IntegratorType integrator;
  FieldPortalConstType field;
};


template <typename IntegratorType, typename FieldType, vtkm::IdComponent Size, typename DeviceAdapterTag>
class ParticleAdvectionWorklet
{
public:
  typedef vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>> FieldHandle;
  typedef typename FieldHandle::template ExecutionTypes<DeviceAdapterTag>::PortalConst
    FieldPortalConstType;
  typedef ParticleAdvectWorklet<IntegratorType,
                                FieldType,
                                Size,
                                DeviceAdapterTag>
    ParticleAdvectWorkletType;

  ParticleAdvectionWorklet(const IntegratorType& it)
    : integrator(it)
  {

  }

  template <typename PointStorage, typename FieldStorage>
  vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>, PointStorage> Run(
    const vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>, PointStorage>& _pl,
      const vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>, PointStorage>& _pr,
    const vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>, FieldStorage> fieldArray)
  {
    pl = _pl;
    pr = _pr;
    field = fieldArray.PrepareForInput(DeviceAdapterTag());
    return run();
  }

  ~ParticleAdvectionWorklet() {}

private:
  vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, 2>> run(bool dumpOutput = false)
  {
    typedef typename vtkm::worklet::DispatcherMapField<ParticleAdvectWorkletType>
      ParticleWorkletDispatchType;


    ParticleAdvectWorkletType particleWorklet(integrator, field);
    ParticleWorkletDispatchType particleWorkletDispatch(particleWorklet);
    particleWorkletDispatch.Invoke(pl, pr);


    return pr;
  }

  IntegratorType integrator;
  vtkm::cont::ArrayHandle<vtkm::Vec<FieldType, Size>> pl, pr;
  vtkm::cont::DataSet ds;
  vtkm::Id maxSteps;
  vtkm::Id ParticlesPerRound;
  FieldPortalConstType field;
};




#endif // vtk_m_worklet_particleadvection_ParticleAdvectionWorklets_h
