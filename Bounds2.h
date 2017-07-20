//============================================================================
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2016 Sandia Corporation.
//  Copyright 2016 UT-Battelle, LLC.
//  Copyright 2016 Los Alamos National Security.
//
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
//  Laboratory (LANL), the U.S. Government retains certain rights in
//  this software.
//============================================================================

#ifndef Bounds_h
#define Bounds_h

#include <vtkm/Range.h>

struct Bounds
{
  vtkm::Range X;
  vtkm::Range Y;

  VTKM_EXEC_CONT
  Bounds() {}

  VTKM_EXEC_CONT
  Bounds(const vtkm::Range& xRange, const vtkm::Range& yRange)
    : X(xRange)
    , Y(yRange)
  {
  }

  template <typename T1, typename T2, typename T3, typename T4>
  VTKM_EXEC_CONT Bounds(const T1& minX,
                        const T2& maxX,
                        const T3& minY,
                        const T4& maxY)
    : X(vtkm::Range(minX, maxX))
    , Y(vtkm::Range(minY, maxY))
  {
  }

  /// Initialize bounds with an array of 6 values in the order xmin, xmax,
  /// ymin, ymax, zmin, zmax.
  ///
  template <typename T>
  VTKM_EXEC_CONT explicit Bounds(const T bounds[4])
    : X(vtkm::Range(bounds[0], bounds[1]))
    , Y(vtkm::Range(bounds[2], bounds[3]))

  {
  }

  /// Initialize bounds with the minimum corner point and the maximum corner
  /// point.
  ///
  template <typename T>
  VTKM_EXEC_CONT Bounds(const vtkm::Vec<T, 2>& minPoint, const vtkm::Vec<T, 2>& maxPoint)
    : X(vtkm::Range(minPoint[0], maxPoint[0]))
    , Y(vtkm::Range(minPoint[1], maxPoint[1]))
  {
  }

  VTKM_EXEC_CONT
  const Bounds& operator=(const Bounds& src)
  {
    this->X = src.X;
    this->Y = src.Y;
    return *this;
  }

  /// \b Determine if the bounds are valid (i.e. has at least one valid point).
  ///
  /// \c IsNonEmpty returns true if the bounds contain some valid points. If
  /// the bounds are any real region, even if a single point or it expands to
  /// infinity, true is returned.
  ///
  VTKM_EXEC_CONT
  bool IsNonEmpty() const
  {
    return (this->X.IsNonEmpty() && this->Y.IsNonEmpty());
  }

  /// \b Determines if a point coordinate is within the bounds.
  ///
  template <typename T>
  VTKM_EXEC_CONT bool Contains(const vtkm::Vec<T, 2>& point) const
  {
    return (this->X.Contains(point[0]) && this->Y.Contains(point[1]));
  }

  /// \b Returns the center of the range.
  ///
  /// \c Center computes the point at the middle of the bounds. If the bounds
  /// are empty, the results are undefined.
  ///
  VTKM_EXEC_CONT
  vtkm::Vec<vtkm::Float64, 2> Center() const
  {
    return vtkm::Vec<vtkm::Float64, 2>(this->X.Center(), this->Y.Center());
  }

  /// \b Expand bounds to include a point.
  ///
  /// This version of \c Include expands the bounds just enough to include the
  /// given point coordinates. If the bounds already include this point, then
  /// nothing is done.
  ///
  template <typename T>
  VTKM_EXEC_CONT void Include(const vtkm::Vec<T, 2>& point)
  {
    this->X.Include(point[0]);
    this->Y.Include(point[1]);
  }

  /// \b Expand bounds to include other bounds.
  ///
  /// This version of \c Include expands these bounds just enough to include
  /// that of another bounds. Esentially it is the union of the two bounds.
  ///
  VTKM_EXEC_CONT
  void Include(const Bounds& bounds)
  {
    this->X.Include(bounds.X);
    this->Y.Include(bounds.Y);
  }

};


/// Helper function for printing bounds during testing
///
static inline VTKM_CONT std::ostream& operator<<(std::ostream& stream, const Bounds& bounds)
{
  return stream << "{ X:" << bounds.X << ", Y:" << bounds.Y << " }";
}

#endif //vtk_m_Bounds_h
