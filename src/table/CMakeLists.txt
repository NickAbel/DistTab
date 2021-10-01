# Lookup table library
SET(KINDP_src ${SRCTYPE}/kind_params.f90)
SET(TABLE_src ${SRCTABLE}/disttab_table.f90)
SET(LOCALPILE_src ${SRCTABLE}/disttab_local_pile.f90)
ADD_LIBRARY(${TABLELIB} STATIC ${TABLE_src} ${KINDP_src} ${LOCALPILE_src})
