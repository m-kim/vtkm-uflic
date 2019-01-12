//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2015 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//  Copyright 2015 UT-Battelle, LLC.
//  Copyright 2015 Los Alamos National Security.
//
//  Under the terms of Contract DE-NA0003525 with NTESS,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================
#ifndef vtk_m_exec_RandomArray_h
#define vtk_m_exec_RandomArray_h

#include <vtkm/ListTag.h>
#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/ExecutionObjectBase.h>
#include "DeviceAdapterAlgorithmGeneral.h"
namespace vtkm
{
namespace exec
{

/// \brief A type list containing types that can be used with an RandomArray.
///
struct RandomArrayTypeListTag : vtkm::ListTagBase<vtkm::Int32, vtkm::Int64>
{
};

/// A class that can be used to Randomally operate on an array of values safely
/// across multiple instances of the same worklet. This is useful when you have
/// an algorithm that needs to accumulate values in parallel, but writing out a
/// value per worklet might be memory prohibitive.
///
/// To construct an RandomArray you will need to pass in an
/// vtkm::cont::ArrayHandle that is used as the underlying storage for the
/// RandomArray
///
/// Supported Operations: add / compare and swap (CAS)
///
/// Supported Types: 32 / 64 bit signed integers
///
///
template <typename T, typename DeviceAdapterTag>
class RandomArray : public vtkm::cont::ExecutionObjectBase
{
public:
  using ValueType = T;

  VTKM_CONT
  RandomArray()
    : RandomImplementation((vtkm::cont::ArrayHandle<T>()))
  {
  }

  template <typename StorageType>
  VTKM_CONT RandomArray(vtkm::cont::ArrayHandle<T, StorageType> handle)
    : RandomImplementation(handle)
  {
  }

  VTKM_EXEC
  T Set(vtkm::Id index, T a, T b) const
  {
    return this->RandomImplementation.Set(index, a, b);
  }

private:
  DeviceAdapterRandomArrayImplementation<T, DeviceAdapterTag> RandomImplementation;
};
}
} // namespace vtkm::exec

#endif //vtk_m_exec_RandomArray_h
