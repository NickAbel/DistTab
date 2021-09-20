# Fortran floating point kind declarations
SET(TYPE_src ${SRCTYPE}/kind_params.f90)
ADD_LIBRARY(${TYPELIB} STATIC ${TYPE_src})
