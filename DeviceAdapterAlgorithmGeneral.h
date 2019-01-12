//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2014 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//  Copyright 2014 UT-Battelle, LLC.
//  Copyright 2014 Los Alamos National Security.
//
//  Under the terms of Contract DE-NA0003525 with NTESS,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================

#ifndef DeviceAdapterAlgorithmGeneral_h
#define DeviceAdapterAlgorithmGeneral_h

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/ArrayHandleDiscard.h>
#include <vtkm/cont/ArrayHandleImplicit.h>
#include <vtkm/cont/ArrayHandleIndex.h>
#include <vtkm/cont/ArrayHandleStreaming.h>
#include <vtkm/cont/ArrayHandleZip.h>
#include <vtkm/cont/internal/FunctorsGeneral.h>

#include <vtkm/exec/internal/ErrorMessageBuffer.h>
#include <vtkm/exec/internal/TaskSingular.h>

#include <vtkm/BinaryPredicates.h>
#include <vtkm/TypeTraits.h>

#include <vtkm/internal/Windows.h>

#include <type_traits>
#include <cstdlib>

#include <vtkm/cont/DeviceAdapterAlgorithm.h>

/// \brief Class providing a device-specific random interface.
///
/// The class provide the actual implementation used by vtkm::exec::AtomicArray.
/// A serial default implementation is provided. But each device will have a different
/// implementation.
///
/// Serial requires no form of atomicity
///
template <typename T, typename DeviceTag>
class DeviceAdapterRandomArrayImplementation
{
public:
  VTKM_CONT
  DeviceAdapterRandomArrayImplementation(
    vtkm::cont::ArrayHandle<T, vtkm::cont::StorageTagBasic> handle)
    : Iterators(IteratorsType(handle.PrepareForInPlace(DeviceTag())))
  {
  }

  inline  T Set(vtkm::Id index, T a, T b) const
  {
      T *lockedValue;
#if defined(_ITERATOR_DEBUG_LEVEL) && _ITERATOR_DEBUG_LEVEL > 0
      using IteratorType = typename vtkm::cont::ArrayPortalToIterators<PortalType>::IteratorType;
      typename IteratorType::pointer temp =
          &(*(Iterators.GetBegin() + static_cast<std::ptrdiff_t>(index)));
      lockedValue = temp;
#else
      lockedValue = (Iterators.GetBegin() + index);
#endif

    return vtkmSet(lockedValue, index, a, b);
  }

private:
  T a,b;
  using PortalType =
    typename vtkm::cont::ArrayHandle<T, vtkm::cont::StorageTagBasic>::template ExecutionTypes<
      DeviceTag>::Portal;
  using IteratorsType = vtkm::cont::ArrayPortalToIterators<PortalType>;
  IteratorsType Iterators;
  inline vtkm::Int64 vtkmSet(vtkm::Int32 * address,
                              const vtkm::Id & idx,
                             T a,
                             T b) const
  {

    address[0] = rand() % b;
    return address[0];
//    return atomicCAS((unsigned long long int*)address,
//                     (unsigned long long int)oldValue,
//                     (unsigned long long int)newValue);
  }
};


#endif //DeviceAdapterAlgorithmGeneral_h
