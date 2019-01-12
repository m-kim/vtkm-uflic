//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2016 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
//  Copyright 2016 UT-Battelle, LLC.
//  Copyright 2016 Los Alamos National Security.
//
//  Under the terms of Contract DE-NA0003525 with NTESS,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================
#ifndef TypeCheckTagRandomArray_h
#define TypeCheckTagRandomArray_h

#include <vtkm/cont/arg/TypeCheck.h>

#include <vtkm/ListTag.h>

#include <vtkm/cont/ArrayHandle.h>
#include <vtkm/cont/StorageBasic.h>

#include "RandomArray.h"

/// The Random array type check passes for an \c ArrayHandle of a structure
/// that is valid for Random access. There are many restrictions on the
/// type of data that can be used for an Random array.
///
template <typename TypeList = vtkm::exec::RandomArrayTypeListTag>
struct TypeCheckTagRandomArray
{
  VTKM_IS_LIST_TAG(TypeList);
};

template <typename TypeList, typename ArrayType>
struct vtkm::cont::arg::TypeCheck<TypeCheckTagRandomArray<TypeList>, ArrayType>
{
  static const bool value = false;
};

template <typename T, typename TypeList>
struct vtkm::cont::arg::TypeCheck<TypeCheckTagRandomArray<TypeList>,
                 vtkm::cont::ArrayHandle<T, vtkm::cont::StorageTagBasic>>
{
  static const bool value = (vtkm::ListContains<TypeList, T>::value &&
                             vtkm::ListContains<vtkm::exec::RandomArrayTypeListTag, T>::value);
};

#endif //TypeCheckTagRandomArray_h
