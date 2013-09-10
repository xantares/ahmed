#set(CMAKE_VERBOSE_MAKEFILE on)
cmake_minimum_required(VERSION 2.6)

# Modules
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)
#include(CheckCXXSourceRuns)

############################################################################
### options

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
   message(STATUS "Setting build type to 'Release' as none was specified.")
   set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
   # Set the possible values of build type for cmake-gui
   set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release")
endif()

project(AHMED)

if(DEFINED ENABLE_64BIT)
   set(ENABLE_64BIT ${ENABLE_64BIT} CACHE BOOL "Enable/disable 64 bit")
else()
   set(ENABLE_64BIT TRUE CACHE BOOL "Enable/disable 64 bit")
endif()

set(BUILD_SHARED_LIBS "ON" CACHE BOOL "Build shared/static libraries")

############################################################################
### compiler flags

# OpenMP
find_package(OpenMP REQUIRED)
if(NOT OPENMP_FOUND)
   message(WARNING "Could not find working OpenMP. Please set manually.")
else()
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Release")
   if(NOT ENABLE_64BIT)
      execute_process(COMMAND ${CMAKE_SOURCE_DIR}/cmake/compilerflags
                              -32 --comp ${CMAKE_C_COMPILER}
                      OUTPUT_VARIABLE OPT_FLAGS_)
   else()
      execute_process(COMMAND ${CMAKE_SOURCE_DIR}/cmake/compilerflags
                              -64 --comp ${CMAKE_C_COMPILER}      
                      OUTPUT_VARIABLE OPT_FLAGS_)
   endif()
else()
   if(NOT ENABLE_64BIT)
      execute_process(COMMAND ${CMAKE_SOURCE_DIR}/cmake/compilerflags
                              --opt -32 --comp ${CMAKE_C_COMPILER}
                      OUTPUT_VARIABLE OPT_FLAGS_)
   else()
      execute_process(COMMAND ${CMAKE_SOURCE_DIR}/cmake/compilerflags
                              --opt -64 --comp ${CMAKE_C_COMPILER}      
                      OUTPUT_VARIABLE OPT_FLAGS_)
  endif()
endif()
string(REGEX REPLACE "\n" "" OPT_FLAGS "${OPT_FLAGS_}")
message(STATUS "OPT_FLAGS " ${OPT_FLAGS})
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OPT_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OPT_FLAGS}")

############################################################################
### external libraries

### METIS
set(ENABLE_METIS "TRUE" CACHE BOOL "Enable/disable METIS")
find_package(Metis REQUIRED)
if(METIS_FOUND)
   set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DMETIS_VERSION=${METIS_VERSION}")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DMETIS_VERSION=${METIS_VERSION}")
endif()

### MPI
set(ENABLE_MPI "" CACHE BOOL "Enable/disable MPI")
if(ENABLE_MPI)
   find_program(MPI_CXX NAMES openmpic++ mpic++ mpicxx mpiCC hcp
                              mpxlC_r mpxlC mpCC cmpic++)
#   message(STATUS "MPI_CXX" ${MPI_CXX}) 
   if(MPI_CXX_FOUND STREQUAL "MPI_CXX-NOTFOUND")
      message(FATAL_ERROR "Could not find MPI C++.")
   endif()
   execute_process(COMMAND ${MPI_CXX} --showme:compile
                   OUTPUT_VARIABLE OUT_VAR)
   string(REGEX REPLACE "-I" "" MPI_CXX_INCLUDE_DIR_ ${OUT_VAR})
   string(REGEX REPLACE "\n" " " MPI_CXX_INCLUDE_DIR "${MPI_CXX_INCLUDE_DIR_}")
   execute_process(COMMAND ${MPI_CXX} --showme:link
                   OUTPUT_VARIABLE MPI_CXX_LIBRARIES_)
   string(REGEX REPLACE "\n" "" MPI_CXX_LIBRARIES
          "${MPI_CXX_LIBRARIES_}")

endif()

### BLAS / LAPACK
set(BLAS_LIBS "-lblas" CACHE STRING "BLAS libraries")
set(LAPACK_LIBS "-llapack" CACHE STRING "LAPACK libraries")

############################################################################
### check whether BLAS works

file(READ "${CMAKE_SOURCE_DIR}/cmake/test_blas.cpp" BLASTEST)
file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_blas/src.cxx" "${BLASTEST}\n")

message(STATUS "Checking whether BLAS works...")

try_run(BLAS_RESULT BLAS_RESULT_COMPILED
        ${CMAKE_BINARY_DIR}
	${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_blas/src.cxx
	CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${BLAS_LIBS}"
        COMPILE_OUTPUT_VARIABLE comp)
#message(STATUS "Output " ${comp})

if(BLAS_RESULT_COMPILED)
   if ("${BLAS_RESULT}" EQUAL 0)
      unset(BLAS_RESULT CACHE)
      unset(BLAS_RESULT_COMPILED CACHE)
      message(STATUS "BLAS compiles and runs.")
   else()
      unset(BLAS_RESULT CACHE)
      unset(BLAS_RESULT_COMPILED CACHE)
      message(FATAL_ERROR "BLAS compiles but test does not execute successfully.")
   endif()
else()
   unset(BLAS_RESULT CACHE)
   unset(BLAS_RESULT_COMPILED CACHE)
   message(FATAL_ERROR "BLAS test does not compile. Check libraries.")
endif()


############################################################################
### check whether LAPACK works

file(READ "${CMAKE_SOURCE_DIR}/cmake/test_lapack.cpp" LAPACKTEST)
file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_lapack/src.cxx" "${LAPACKTEST}\n")

set(LAPACK_BLAS_LIBS "${LAPACK_LIBS};${BLAS_LIBS}")

message(STATUS "Checking whether LAPACK works...")

try_run(LAPACK_RESULT LAPACK_RESULT_COMPILED ${CMAKE_BINARY_DIR}
	${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_lapack/src.cxx
	CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${LAPACK_BLAS_LIBS}"
	RUN_OUTPUT_VARIABLE LAPACK_OUTPUT)
if(LAPACK_RESULT_COMPILED)
   if("${LAPACK_RESULT}" EQUAL 0)
      unset(LAPACK_RESULT CACHE)
      unset(LAPACK_RESULT_COMPILED CACHE)
#      message(STATUS ${LAPACK_OUTPUT})
      string(REGEX REPLACE "[a-zA-Z\n ]" "" LAPACK_VERSION "${LAPACK_OUTPUT}")
      message(STATUS "LAPACK compiles and runs.")
   else()
      unset(LAPACK_RESULT CACHE)
      unset(LAPACK_RESULT_COMPILED CACHE)
      message(FATAL_ERROR "LAPACK compiles but test does not execute successfully")
   endif()
else()
   unset(LAPACK_RESULT CACHE)
   unset(LAPACK_RESULT_COMPILED CACHE)
   message(FATAL_ERROR "LAPACK test does not compile. Check libraries")
endif()


############################################################################
### find correct calling convention for dotc

file(READ "${CMAKE_SOURCE_DIR}/cmake/test_zdotc.cpp" ZDOTCTEST)
file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_zdotc/src.cxx" "${ZDOTCTEST}\n")

message(STATUS "Test if ZDOTC works with GFORTRAN calling convention...")

try_run(ZDOT_RESULT ZDOT_RESULT_COMPILED
        ${CMAKE_BINARY_DIR}
	${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_zdotc/src.cxx
	CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${BLAS_LIBS}"
	RUN_OUTPUT_VARIABLE ZDOT_OUTPUT)

if(ZDOT_RESULT_COMPILED)
   message(STATUS "Test Output: ${ZDOT_OUTPUT}")
   if ("${ZDOT_RESULT}" EQUAL 0)
      message(STATUS "GFORTRAN calling convention detected")
   else()
      message(STATUS "GFORTRAN calling convention failed. Test for G77 calling convention")
      file(WRITE "${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_zdotc/src.cxx" "${ZDOTCTEST}\n")
      try_run(ZDOTG77_RESULT ZDOTG77_RESULT_COMPILED ${CMAKE_BINARY_DIR}
	      ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeTmp/test_zdotc/src.cxx
              CMAKE_FLAGS "-DLINK_LIBRARIES:STRING=${BLAS_LIBS}"
	      COMPILE_DEFINITIONS "-DG77CONVENTION"
	      RUN_OUTPUT_VARIABLE ZDOT_OUTPUT)
      message(STATUS "Test Output: ${ZDOT_OUTPUT}")
      if ("${ZDOTG77_RESULT}" EQUAL 0)
         message(STATUS "G77 calling convention detected")
	 set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DG77CONVENTION")
      else()
         message(WARNING "AHMED will likely not work in complex mode.")
      endif()
   endif()
else()
   unset(ZDOT_RESULT CACHE)
   unset(ZDOT_RESULT_COMPILED CACHE)
   message(FATAL_ERROR "ZDOTC test does not compile. Check libraries")
endif()

unset(ZDOT_RESULT CACHE)
unset(ZDOT_RESULT_COMPILED CACHE)
unset(ZDOTG77_RESULT CACHE)
unset(ZDOTG77_RESULT_COMPILED CACHE)

############################################################################
### report summary

message(STATUS
        "----------------------------------------------------------------")
message(STATUS "                               Summary")
message(STATUS
        "----------------------------------------------------------------")

message(STATUS "C compiler                 " ${CMAKE_C_COMPILER})
message(STATUS "C compiler flags           " ${CMAKE_C_FLAGS})
message(STATUS "C++ compiler               " ${CMAKE_CXX_COMPILER})
message(STATUS "C++ compiler flags         " ${CMAKE_CXX_FLAGS})
message(STATUS "BLAS library               " ${BLAS_LIBS})
message(STATUS "LAPACK library             " ${LAPACK_LIBS})
message(STATUS "LAPACK version             " ${LAPACK_VERSION})

if(ENABLE_METIS)
   message(STATUS "METIS include directory    " ${METIS_INCLUDE_DIR})
   message(STATUS "METIS library              " ${METIS_LIBRARY})
   message(STATUS "METIS version              " ${METIS_VERSION})
else()
   message(STATUS "METIS disabled.")
endif()

if(ENABLE_MPI)
   message(STATUS "MPI C++ compiler           " ${MPI_CXX})
   message(STATUS "MPI include directory      " ${MPI_CXX_INCLUDE_DIR})
   message(STATUS "MPI library path           " ${MPI_CXX_LIBRARIES})
else()
   message(STATUS "MPI disabled.")
endif()

message(STATUS "Build type                 " ${CMAKE_BUILD_TYPE})

if(ENABLE_64BIT)
   message(STATUS "Architecture               64 bit" )
else()
   message(STATUS "Architecture               32 bit" )
endif()

if(BUILD_SHARED_LIBS)
   message(STATUS "Building shared libraries.")
else()
   message(STATUS "Building static libraries.")
endif()
message(STATUS
        "----------------------------------------------------------------")


############################################################################
### add files to library

file(GLOB H_CPP H/HLU.cpp H/mltaGeHHeH.cpp H/Hequilib.cpp H/mltaGeHLtH.cpp
                H/Htrunc.cpp H/mltaGeHUtHh.cpp H/TH_solve.cpp H/mltaGeHVec.cpp
                H/addHH.cpp H/mltaGeHVec_omp.cpp H/agglH.cpp H/mblock.cpp
                H/mltaGeHhDiHGeH.cpp H/convCS_toH.cpp H/mltaUtHGeH.cpp
                H/mltaGeHhDiHGeH_toHeH.cpp H/initH.cpp H/mltaGeHhGeH.cpp
                H/mltaGeHhGeH_toHeH.cpp H/mblock_D.cpp H/mblock_S.cpp
                H/mltaHeHGeH.cpp H/mblock_C.cpp H/mltaLtHGeH.cpp
                H/mltaUtHUtHh.cpp H/mblock_Z.cpp H/mltaUtHhGeH.cpp
                H/mltaGeHGeH.cpp H/mltaUtHhUtH_toHeH.cpp H/mltaGeHGeHh.cpp H/nrmH.cpp
                H/mltaGeHGeHh_toHeH.cpp H/psoutH.cpp)

file(GLOB BASMOD_CPP basmod/progress.cpp)

file(GLOB BASMOD_C basmod/cputime.c basmod/realtime.c)

file(GLOB SOLVERS_CPP solvers/BiCGStab.cpp solvers/CG.cpp
                      solvers/GMRES.cpp solvers/FGMRES.cpp solvers/MINRES.cpp)

file(GLOB SPARSE_CPP sparse/CS_CRS2CRSSym.cpp sparse/CS_perm.cpp
                     sparse/CS_CRSSym2CRS.cpp sparse/CS_gen.cpp
                     sparse/CS_genDiagPrec.cpp sparse/CS_getDist.cpp
                     sparse/amux.cpp sparse/amux_symm.cpp
                     sparse/atmux.cpp sparse/CS_getRemotlyNode.cpp
                     sparse/CS_io.cpp sparse/CS_is_connected.cpp
                     sparse/genAdjM.cpp sparse/CS_transp.cpp)

file(GLOB CLUSTER_CPP cluster/AdjMatrix.cpp cluster/binSearch.cpp
                      cluster/ClusterAlg.cpp cluster/specialSort.cpp
                      cluster/CS_getPathAndDist.cpp cluster/blcluster.cpp)
if(ENABLE_METIS)
   list(APPEND CLUSTER_CPP "${CMAKE_SOURCE_DIR}/cluster/cluster_alg.cpp" 
   	                   "${CMAKE_SOURCE_DIR}/cluster/Separator.cpp")
endif()

file(GLOB PARALLEL_CPP parallel/gen_blseq.cpp parallel/parallel.cpp)
if(ENABLE_MPI)
   list(APPEND PARALLEL_CPP "${CMAKE_SOURCE_DIR}/parallel/CG_MPI.cpp" 
                            "${CMAKE_SOURCE_DIR}/parallel/GMRES_MPI.cpp"
                            "${CMAKE_SOURCE_DIR}/parallel/psoutH_MPI.cpp"
                            "${CMAKE_SOURCE_DIR}/parallel/mltaGeHVec_MPI.cpp"
                            "${CMAKE_SOURCE_DIR}/parallel/transfH_MPI.cpp")

   file(GLOB ND_CPP ND/Blas_ND.cpp ND/HUhDU.cpp ND/mltaGeHVec_ND.cpp
                    ND/CG_ND.cpp ND/Helpers_ND.cpp ND/TH_solve.cpp
                    ND/multaCRSvec_ND.cpp ND/HCholesky_ND.cpp
                    ND/subdivCRS_ND.cpp HD/HLU_ND.cpp ND/convCS_toH_ND.cpp)
endif()

file(GLOB LIB_INCLUDES Include/*.h)

# add includes
include_directories("${CMAKE_SOURCE_DIR}/Include")

if(ENABLE_METIS)
   include_directories("${METIS_INCLUDE_DIR}")
endif()

if(ENABLE_MPI)
   include_directories("${MPI_CXX_INCLUDE_DIR}")
endif()

# add library target
add_library(AHMED ${H_CPP} ${BASMOD_CPP} ${BASMOD_C} ${SOLVERS_CPP} ${CLUSTER_CPP}
                  ${SPARSE_CPP} ${PARALLEL_CPP} ${ND_CPP} ${LIB_INCLUDES})

############################################################################
### link library

target_link_libraries(AHMED ${BLAS_LIBS} ${LAPACK_LIBS})

if(ENABLE_METIS)
   target_link_libraries(AHMED ${METIS_LIBRARY})
endif()

if(ENABLE_MPI)
   target_link_libraries(AHMED ${MPI_CXX_LIBRARIES})
endif()

############################################################################
### install library

if (WIN32)
  set (INSTTARGET "$ENV{ProgramFiles}/ahmed/")
else (WIN32)
  set (INSTTARGET "/usr/local")
endif (WIN32)

install(TARGETS AHMED
	RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
     
# Install headers

if (NOT ENABLE_METIS)
   list(REMOVE_ITEM LIB_INCLUDES "${CMAKE_SOURCE_DIR}/Include/cluster_alg.h")
endif ()

install(FILES ${LIB_INCLUDES} DESTINATION ${CMAKE_INSTALL_PREFIX}/include/AHMED)
