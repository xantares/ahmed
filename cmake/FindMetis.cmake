#
# Find METIS headers and lib
#
# METIS_INCLUDE_DIR - where to find metis.h
# METIS_LIBRARY     - libmetis.a

if(DEFINED METIS_ROOT)
   set(METIS_ROOT ${METIS_ROOT} CACHE BOOL "METIS root path.")
else()
   #check if METIS_ROOT is defined in envir and use it
   if (NOT $ENV{METIS_ROOT} STREQUAL "")
      set(METIS_ROOT $ENV{METIS_ROOT} CACHE STRING "METIS root path.")
   endif()
endif()


IF(ENABLE_METIS)

   # figure out what the next statement does
   IF(METIS_ROOT)
      file(TO_CMAKE_PATH ${METIS_ROOT} METIS_ROOT)
   ENDIF()

   FIND_PATH(METIS_INCLUDE_DIR metis.h
      ${METIS_ROOT}
      ${METIS_ROOT}/include
      ${METIS_ROOT}/Include
      /opt/metis/include
      /usr/local/include
      /usr/include
      /usr/include/metis
   )

   FIND_LIBRARY(METIS_LIBRARY metis
      ${METIS_ROOT}
      ${METIS_ROOT}/Lib
      ${METIS_ROOT}/lib
      /opt/metis/lib
      /usr/local/lib
      /usr/lib
   )


   IF(METIS_INCLUDE_DIR AND METIS_LIBRARY)
      SET(METIS_FOUND "YES")
      # find out version
      execute_process(
        WORKING_DIRECTORY ${METIS_INCLUDE_DIR}
        COMMAND grep METIS_VER_MAJOR metis.h 
        OUTPUT_VARIABLE METIS_VERSION_
      )
      STRING(REGEX REPLACE "[a-zA-Z_ \#\n]" "" METIS_VERSION "${METIS_VERSION_}")
      IF(METIS_VERSION STREQUAL "")
         SET(METIS_VERSION "4")
      ENDIF()
   ENDIF()

   #include(FindPackageHandleStandardArgs)
   #find_package_handle_standard_args(METIS DEFAULT_MSG METIS_INCLUDE_DIR METIS_LIBRARY)

   #MARK_AS_ADVANCED( METIS_INCLUDE_DIR METIS_LIBRARY )
ENDIF()
