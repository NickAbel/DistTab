# Parallel communication
SET(KINDP_src ${SRCTYPE}/kind_params.f90)
SET(PARALLEL_src ${SRCPARALLEL}/communicator.f90)
ADD_LIBRARY(${PARALLELLIB} STATIC ${PARALLEL_src} ${KINDP_src})
