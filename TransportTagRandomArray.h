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
#ifndef TransportTagRandomArray_h
#define TransportTagRandomArray_h

#include <vtkm/Types.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/StorageBasic.h>

#include <vtkm/cont/arg/Transport.h>

#include "RandomArray.h"


/// \brief \c Transport tag for in-place arrays with Random operations.
///
/// \c TransportTagRandomArray is a tag used with the \c Transport class to
/// transport \c ArrayHandle objects for data that is both input and output
/// (that is, in place modification of array data). The array will be wrapped
/// in a vtkm::exec::RandomArray class that provides Random operations (like
/// add and compare/swap).
///
struct TransportTagRandomArray
{
};

template <typename T, typename Device>
struct vtkm::cont::arg::Transport<TransportTagRandomArray,
                 vtkm::cont::ArrayHandle<T, vtkm::cont::StorageTagBasic>,
                 Device>
{
  using ExecObjectType = vtkm::exec::RandomArray<T, Device>;

  template <typename InputDomainType>
  VTKM_CONT ExecObjectType operator()(vtkm::cont::ArrayHandle<T, vtkm::cont::StorageTagBasic> array,
                                      const InputDomainType&,
                                      vtkm::Id,
                                      vtkm::Id) const
  {
    // Note: we ignore the size of the domain because the randomly accessed
    // array might not have the same size depending on how the user is using
    // the array.

    return ExecObjectType(array);
  }
};

#endif //TransportTagRandomArray_h
