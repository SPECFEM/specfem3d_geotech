# make utilities
# REVISION:
#   HNG, Sep 06,2011;

# build exodus2sem
IF(BUILD_UTILITIES_EXODUS2SEM)
  ADD_EXECUTABLE(exodus2sem
              exodus2sem.c)
ENDIF(BUILD_UTILITIES_EXODUS2SEM)

# build write_sos
IF(BUILD_UTILITIES_WRITE_SOS)
  ADD_EXECUTABLE(write_sos
              write_sos.f90)
ENDIF(BUILD_UTILITIES_WRITE_SOS)

