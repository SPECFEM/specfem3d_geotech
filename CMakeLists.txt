# cmake file for the cross-platform configure of the
# software package SPECFEM3D_GEOTECH
# REVISION:
#   HNG, Sep 06,2011;

# required cmake version
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.4 FATAL_ERROR)

# define project
PROJECT(semgeotech Fortran C)

# Disallow in-source build
IF("${semgeotech_SOURCE_DIR}"  STREQUAL "${semgeotech_BINARY_DIR}")
  MESSAGE(FATAL_ERROR "SPECFEM3D_GEOTECH requires an out-of-source build. Please create a separate binary directory and run CMake there.")
ENDIF("${semgeotech_SOURCE_DIR}"  STREQUAL "${semgeotech_BINARY_DIR}")

# Test Fortran compiler
TRY_COMPILE(TEST_F90 ${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/src/testf90.f90)
IF(NOT TEST_F90)
  MESSAGE(FATAL_ERROR "Fortran compiler: ${CMAKE_Fortran_COMPILER} doesn't work! Please set working CMAKE_Fortran_COMPILER!")
ENDIF(NOT TEST_F90)

# build directories
SET(BDIR ${CMAKE_BINARY_DIR}/bin)
SET(DDIR ${CMAKE_BINARY_DIR}/doc)
SET(IDIR ${CMAKE_BINARY_DIR}/input)
SET(LDIR ${CMAKE_BINARY_DIR}/lib)
SET(ODIR ${CMAKE_BINARY_DIR}/output)
SET(PDIR ${CMAKE_BINARY_DIR}/partition)
SET(SDIR ${CMAKE_SOURCE_DIR})
FILE(MAKE_DIRECTORY ${BDIR})
#FILE(MAKE_DIRECTORY ${DDIR})
FILE(MAKE_DIRECTORY ${IDIR}) 
FILE(MAKE_DIRECTORY ${ODIR}) 
FILE(MAKE_DIRECTORY ${PDIR}) 

# set paths to executable and libraries
SET(EXECUTABLE_OUTPUT_PATH ${BDIR})
SET(LIBRARY_OUTPUT_PATH ${LDIR})

# set internal default logical variables
SET(ISDOCUMENTATION OFF)
SET(ISPARTMESH OFF)
SET(ISPARALLEL OFF)

# build documentation option
SET(BUILD_DOCUMENTATION OFF CACHE BOOL "Build documentation for
SPECFEM3D_GEOTECH")
IF(BUILD_DOCUMENTATION)
  SET(ISDOCUMENTATION ON)
  FIND_PACKAGE(LATEX REQUIRED)
  IF(NOT PDFLATEX_COMPILER OR NOT BIBTEX_COMPILER)
    MESSAGE(WARNING "LaTeX compiler/s not found! DOCUMENTATION cannot be built!"    )
  SET(ISDOCUMENTATION OFF)
  ENDIF(NOT PDFLATEX_COMPILER OR NOT BIBTEX_COMPILER)
ENDIF(BUILD_DOCUMENTATION)

#SET(BUILD_SLOPE ON CACHE BOOL "Build slope stability program")
#SET(BUILD_EXCAVATION ON CACHE BOOL "Build excavation program")

# build partmesh option
SET(BUILD_PARTMESH OFF CACHE BOOL "Build mesh partitioning program")
IF(BUILD_PARTMESH)
  SET(ISPARTMESH ON)
  FIND_PATH(SCOTCH_LIBRARY_PATH NAMES libscotch.a libscotcherr.a
  PATHS /usr/lib /usr/local/lib ENV LD_LIBRARY_PATH)
  IF(NOT SCOTCH_LIBRARY_PATH)
    MESSAGE(WARNING "SCOTCH library not found! Please set correct SCOTCH_LIBRARY_PATH manually!")
    SET(ISPARTMESH OFF)
  ENDIF(NOT SCOTCH_LIBRARY_PATH)
ENDIF(BUILD_PARTMESH)
SET(BUILD_UTILITIES_EXODUS2SEM OFF CACHE BOOL "Build exodus2sem program")
SET(BUILD_UTILITIES_WRITE_SOS OFF CACHE BOOL "Build write_sos program")

# MPI for parallel processing
SET(ENABLE_MPI OFF CACHE BOOL "Build semgeotech with MPI support")
IF(ENABLE_MPI)
  SET(ISPARALLEL ON)
  FIND_PACKAGE(MPI COMPONENTS Fortran REQUIRED)
  IF(NOT MPI_FOUND)
    MESSAGE(WARNING "MPI Fortran compiler not found! Parallel semgeotech cannot be built!")
    SET(ISPARALLEL OFF)
  ENDIF(NOT MPI_FOUND)
  
  # Test MPI Fortran compiler
  #include_directories(${MPI_Fortran_INCLUDE_PATH})
  TRY_COMPILE(TEST_MPI ${CMAKE_BINARY_DIR}/CMakeTmp ${CMAKE_SOURCE_DIR}/src/testmpi.f90
  COMPILE_DEFINITIONS ${MPI_Fortran_LINK_FLAGS}
  CMAKE_FLAGS "-DINCLUDE_DIRECTORIES=${MPI_Fortran_INCLUDE_PATH}"
  "-DLINK_LIBRARIES=${MPI_Fortran_LIBRARIES}"
  OUTPUT_VARIABLE test)
  IF(NOT TEST_MPI)
    MESSAGE(WARNING "MPI Fortran compiler: ${CMAKE_Fortran_COMPILER} doesn't work! Parallel semgeotech cannot be built!")
    SET(ISPARALLEL OFF)
  ENDIF(NOT TEST_MPI) 
INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH}) 
ENDIF(ENABLE_MPI)

# process CMakeLists.txt in other subdirectories
ADD_SUBDIRECTORY(doc)
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(utilities)

# notes
#SET_SOURCE_FILES_PROPERTIES(${SRC} PROPERTIES LANGUAGE "C")

