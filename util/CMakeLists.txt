# make utilities
# AUTHOR
#   Hom Nath Gharti
# REVISION:
#   HNG, Sep 06,2011;

# build exodus2semgeotech
IF(BUILD_UTILITIES_EXODUS2SEM)
  ADD_EXECUTABLE(exodus2semgeotech
              exodus2semgeotech.c)
ENDIF(BUILD_UTILITIES_EXODUS2SEM)

# build gid2semgeotech
IF(BUILD_UTILITIES_GID2SEM)
  ADD_EXECUTABLE(gid2semgeotech
              gid2semgeotech.c)
ENDIF(BUILD_UTILITIES_GID2SEM)

# build write_sos
IF(BUILD_UTILITIES_WRITE_SOS)
  ADD_EXECUTABLE(write_sos
              write_sos.f90)
ENDIF(BUILD_UTILITIES_WRITE_SOS)

