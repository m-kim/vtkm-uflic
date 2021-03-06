##=============================================================================
##
##  Copyright (c) Kitware, Inc.
##  All rights reserved.
##  See LICENSE.txt for details.
##
##  This software is distributed WITHOUT ANY WARRANTY; without even
##  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
##  PURPOSE.  See the above copyright notice for more information.
##
##  Copyright 2015 Sandia Corporation.
##  Copyright 2015 UT-Battelle, LLC.
##  Copyright 2015 Los Alamos National Security.
##
##  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
##  the U.S. Government retains certain rights in this software.
##  Under the terms of Contract DE-AC52-06NA25396 with Los Alamos National
##  Laboratory (LANL), the U.S. Government retains certain rights in
##  this software.
##
##=============================================================================
cmake_minimum_required(VERSION 3.3)
project(vtkm-uflic)


#Find the VTK-m package
find_package(VTKm REQUIRED
  OPTIONAL_COMPONENTS OpenMP Serial TBB CUDA
  )


set(UFLIC_SRC
  screenspace.cpp
)
set(UFLIC_HDR
    DeviceAdapterAlgorithmCuda.h
    DeviceAdapterAlgorithmGeneral.h
    RandomArray.h
    TransportTagRandomArray.h
    TypeCheckTagRandomArray.h
    RandomMapField.h
  Bounds2.h
  Convolve.h
  Draw.h
  Evaluator.h
  Integrator.h
  Jitter.h
  Normalize.h
  Reader.h
  Sharpen.h
  ParticleAdvection.h
)


add_executable(UFLIC_rendering_Serial ${UFLIC_SRC} ${UFLIC_HDR})
#target_include_directories(UFLIC_rendering_Serial PRIVATE vtkm_cont )
target_link_libraries(UFLIC_rendering_Serial vtkm_cont vtkm_rendering )
target_compile_options(UFLIC_rendering_Serial PRIVATE ${VTKm_COMPILE_OPTIONS})
target_compile_definitions(UFLIC_rendering_Serial PRIVATE
                       "VTKM_DEVICE_ADAPTER=VTKM_DEVICE_ADAPTER_SERIAL")

if(VTKm_ENABLE_CUDA)
  set(old_nvcc_flags ${CUDA_NVCC_FLAGS})
  set(old_cxx_flags ${CMAKE_CXX_FLAGS})
  vtkm_setup_nvcc_flags( old_nvcc_flags old_cxx_flags)
  vtkm_disable_troublesome_thrust_warnings()

  # Cuda compiles do not respect target_include_directories
  cuda_include_directories(${VTKm_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR})

  set (cudaSource  "${CMAKE_CURRENT_BINARY_DIR}/main.cu")
  configure_file(main.cpp "${CMAKE_CURRENT_BINARY_DIR}/main.cu" COPYONLY)
  cuda_add_executable(UFLIC_rendering_CUDA ${cudaSource})
  target_include_directories(UFLIC_rendering_CUDA PRIVATE vtkm_cont ${CMAKE_CURRENT_SOURCE_DIR})
  target_link_libraries(UFLIC_rendering_CUDA PRIVATE vtkm::cuda)
  target_compile_options(UFLIC_rendering_CUDA PRIVATE ${VTKm_COMPILE_OPTIONS})
   set(CUDA_NVCC_FLAGS ${old_nvcc_flags})
  set(CMAKE_CXX_FLAGS ${old_cxx_flags})
endif()

if(VTKm_ENABLE_TBB)
    add_executable(UFLIC_rendering_TBB ${UFLIC_SRC} ${UFLIC_HDR})
    target_include_directories(UFLIC_rendering_TBB PRIVATE vtkm_cont)
    target_link_libraries(UFLIC_rendering_TBB PRIVATE vtkm::tbb)
    target_compile_options(UFLIC_rendering_TBB PRIVATE ${VTKm_COMPILE_OPTIONS})
    target_compile_definitions(UFLIC_rendering_TBB PRIVATE "VTKM_DEVICE_ADAPTER=VTKM_DEVICE_ADAPTER_TBB")

endif()

if(VTKm_ENABLE_OPENMP)
    add_executable(UFLIC_rendering_OpenMP ${UFLIC_SRC} ${UFLIC_HDR})
    #    target_include_directories(UFLIC_rendering_OpenMP PRIVATE vtkm_cont)
    target_link_libraries(UFLIC_rendering_OpenMP PRIVATE vtkm_cont vtkm_rendering)
    target_compile_options(UFLIC_rendering_OpenMP PRIVATE ${VTKm_COMPILE_OPTIONS})
    target_compile_definitions(UFLIC_rendering_OpenMP PRIVATE "VTKM_DEVICE_ADAPTER=VTKM_DEVICE_ADAPTER_OPENMP")

endif()


