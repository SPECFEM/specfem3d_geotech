# build main programs
# AUTHOR
#   Hom Nath Gharti
# REVISION:
#   HNG, Sep 06,2011;

# build sem library
SET(srcsem
   global.f90
   string_library.f90
   math_library.f90
   shape_library.f90
   element.f90
   integration.f90
   elastic.f90
   read_input.f90
   gll_library.f90
   mesh_spec.f90
   apply_bc.f90
   preprocess.f90
   apply_traction.f90
   compute_pressure.f90
   plastic_library.f90
   visual.f90
   postprocess.f90
   excavation_library.f90)
#IF(BUILD_EXCAVATION)
#  LIST(APPEND srcsem excavation_library.f90)
#ENDIF(BUILD_EXCAVATION)
ADD_LIBRARY(sem ${srcsem})

# build parallel library
IF(ISPARALLEL)
  SET(srcparallel
     set_precision_mpi.f90
     mpi_library.f90
     math_library_mpi.f90
     ghost_library_mpi.f90
     solver_mpi.f90)
  #SET_PROPERTY(SOURCE ${srcparallel} PROPERTY COMPILE_DEFINITIONS USE_MPI)
  ADD_LIBRARY(parallel ${srcparallel})
  ADD_DEPENDENCIES(parallel sem)
  #TARGET_LINK_LIBRARIES(parallel ${MPI_Fortran_LIBRARIES})
ELSE(ISPARALLEL)
  SET(srcserial
     serial_library.f90
     math_library_serial.f90
     solver.f90)
  ADD_LIBRARY(serial ${srcserial})
  ADD_DEPENDENCIES(serial sem)
ENDIF(ISPARALLEL)

# build geotech library
#SET(srcgeotech "")
#message("${srcgeotech}")
#IF(BUILD_EXCAVATION)
#  LIST(APPEND srcgeotech excavation.F90)
#ENDIF(BUILD_EXCAVATION)
#IF(BUILD_SLOPE)
#  LIST(APPEND srcgeotech slope.F90)
#ENDIF(BUILD_SLOPE)
SET(srcgeotech
   excavation.F90
   slope.F90)
#compiler preprocessor directives
#ADD_DEFINITIONS(CMAKE_Fortran_FLAGS "-DUSE_MPI")
IF(ISPARALLEL)
  SET_PROPERTY(SOURCE ${srcgeotech} PROPERTY COMPILE_DEFINITIONS USE_MPI)
  #SET_PROPERTY(SOURCE ${srcgeotech} PROPERTY INCLUDE_DIRECTORIES 
  #                    ${MPI_Fortran_INCLUDE_PATH})
ENDIF(ISPARALLEL)
ADD_LIBRARY(geotech ${srcgeotech})
IF(ISPARALLEL)
  ADD_DEPENDENCIES(geotech sem parallel)
ELSE(ISPARALLEL)
  ADD_DEPENDENCIES(geotech sem serial)
ENDIF(ISPARALLEL)

# build semgeotech
#TARGET_LINK_LIBRARIES(geotech ${semall})
#SET_TARGET_PROPERTIES(geotech PROPERTIES LINK_FLAGS 
#                     "${MPI_Fortran_LINK_FLAGS}")
IF(ISPARALLEL)
  SET_PROPERTY(SOURCE semgeotech.F90 PROPERTY COMPILE_DEFINITIONS USE_MPI)
  #SET_PROPERTY(SOURCE semgeotech.F90 PROPERTY INCLUDE_DIRECTORIES
  #             ${MPI_Fortran_INCLUDE_PATH})
  SET(LIBS geotech parallel sem ${MPI_Fortran_LIBRARIES})
  SET(mainprog psemgeotech)
ELSE(ISPARALLEL)
  SET(LIBS geotech serial sem)
  SET(mainprog semgeotech)
ENDIF(ISPARALLEL)
ADD_EXECUTABLE(${mainprog}
              semgeotech.F90)
TARGET_LINK_LIBRARIES(${mainprog} ${LIBS})
#IF(ISPARALLEL)
#  SET(semall geotech parallel sem)
#ELSE(ISPARALLEL)
#  SET(semall geotech serial sem)
#ENDIF(ISPARALLEL)
#TARGET_LINK_LIBRARIES(semgeotech ${semall} ${MPI_Fortran_LIBRARIES})

# Build partmesh program
# partmesh program partitions the mesh for the parallel processing
IF(ISPARTMESH)
  SET(srcpartmesh
     global.f90
     string_library.f90
     serial_library.f90
     math_library.f90
     partmesh_library.f90
     read_input.f90
     partmesh_scotch.f90
     partmesh.f90
  )
  ADD_LIBRARY(partscotch ${srcpartmesh})
  TARGET_INCLUDE_DIRECTORIES(partscotch PUBLIC
  ${SCOTCH_INCLUDE_PATH})
  ADD_EXECUTABLE(partmesh
              partmesh.f90)
  TARGET_LINK_LIBRARIES(partmesh partscotch
${SCOTCH_LIBRARY_PATH}/libscotch.a ${SCOTCH_LIBRARY_PATH}/libscotcherr.a)
ENDIF(ISPARTMESH)


