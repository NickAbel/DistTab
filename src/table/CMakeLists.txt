# Lookup table library
SET(TABLE_src ${SRCTABLE}/disttab_table.f90)
ADD_LIBRARY(${TABLELIB} STATIC ${TABLE_src})
