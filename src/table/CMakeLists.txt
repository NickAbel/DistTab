# Lookup table library
SET(KINDP_src ${SRCTYPE}/kind_params.f90)
SET(TABLE_src ${SRCTABLE}/disttab_table.f90)
SET(PARALLEL_src ${SRCPARALLEL}/communicator.f90)
ADD_LIBRARY(${TABLELIB} STATIC ${TABLE_src} ${KINDP_src} ${PARALLEL_src})
