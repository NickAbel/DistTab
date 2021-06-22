########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(DISTTAB_src ${SRCTESTS}/disttab_test_partitioning.f90
                ${SRCTESTS}/test.f90 
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${DISTTABEXE} ${DISTTAB_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# This links foo to the bar library
TARGET_LINK_LIBRARIES(${DISTTABEXE} ${TABLELIB})

# Uncomment if you have parallelization
#IF(USE_OPENMP)
#    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
#                          COMPILE_FLAGS "${OpenMP_Fortran_FLAGS}"
#                          LINK_FLAGS "${OpenMP_Fortran_FLAGS}")
#ELSEIF(USE_MPI)
#    SET_TARGET_PROPERTIES(${FOOEXE} PROPERTIES
#                          COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
#                          LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
#    INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
#    TARGET_LINK_LIBRARIES(${FOOEXE} ${MPI_Fortran_LIBRARIES})
#ENDIF(USE_OPENMP)

#####################################
# Tell how to install this executable
#####################################

SET(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
INSTALL(TARGETS ${DISTTABEXE} RUNTIME DESTINATION bin)
