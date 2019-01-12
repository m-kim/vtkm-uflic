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
#ifndef WorkletBase_h
#define WorkletBase_h

#include <vtkm/worklet/WorkletMapField.h>
#include "TypeCheckTagRandomArray.h"
#include "TransportTagRandomArray.h"

class RandomWorklet : public vtkm::worklet::WorkletMapField
{
public:

  /// \c ControlSignature tag for whole input/output arrays.
  ///
  /// The \c RandomArrayInOut control signature tag specifies an \c ArrayHandle
  /// passed to the \c Invoke operation of the dispatcher. This is converted to
  /// a \c vtkm::exec::RandomArray object and passed to the appropriate worklet
  /// operator argument with one of the default args. The provided random
  /// operations can be used to resolve concurrency hazards, but have the
  /// potential to slow the program quite a bit.
  ///
  /// The template operator specifies all the potential value types of the
  /// array. The default value type is all types.
  ///
  template <typename TypeList = AllTypes>
  struct RandomArrayInOut : vtkm::cont::arg::ControlSignatureTagBase
  {
    typedef TypeCheckTagRandomArray<TypeList> TypeCheckTag;
    typedef TransportTagRandomArray TransportTag;
    typedef vtkm::exec::arg::FetchTagExecObject FetchTag;
  };

};
#endif //WorkletBase_h
