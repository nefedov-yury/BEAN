#------------------------------------------------------------------------------
# Identifier: <=> boss/DetectorDescription/Identifier
#------------------------------------------------------------------------------

SET( libname Identifier )

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )

#------------------------------------------------------------------------------
# the list of source files
FILE( GLOB CXX_FILES src/*.cxx )

#------------------------------------------------------------------------------
# do not print warnings about known features
#  IF ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
#       CMAKE_CXX_COMPILER_ID MATCHES "GNU"       )
#    ADD_COMPILE_OPTIONS(-Wsign-compare)
#  ENDIF()

#------------------------------------------------------------------------------
#
IF( CMAKE_SYSTEM_NAME MATCHES Windows ) # WINDOWS
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
ELSE()
  ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )
  MAC_LIB_RPATH( ${libname} )
ENDIF()

#------------------------------------------------------------------------------
# install
#------------------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

