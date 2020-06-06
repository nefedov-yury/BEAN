#------------------------------------------------------------------------------
# AbsCor: <=> boss/Analysis/PhotonCor/AbsCor
#------------------------------------------------------------------------------

SET( libname AbsCor )

INCLUDE_DIRECTORIES(
   ${CMAKE_CURRENT_SOURCE_DIR}
   ${ROOTEVENTDATA_DIR}
   ${BeanCore_DIR}
   ${CLHEP_INCLUDE_DIR}
   ${ROOT_INCLUDE_DIR}
)
IF( BOSS_VERSION  VERSION_GREATER  "7.0.0" )
  INCLUDE_DIRECTORIES( ${Identifier_DIR} )
ENDIF()

#------------------------------------------------------------------------------
# the list of source files
FILE( GLOB CXX_FILES src/*.cxx )

#------------------------------------------------------------------------------
#
IF( CMAKE_SYSTEM_NAME MATCHES Windows ) # WINDOWS
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
ELSE()
  ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )
  MAC_LIB_RPATH( ${libname} )
ENDIF()

TARGET_LINK_LIBRARIES( ${libname}
   RootEventData
   BeanCore
   ${ROOT_LIBRARIES}
   ${CLHEP_LIBRARIES}
)

IF( BOSS_VERSION  VERSION_GREATER  "7.0.0" )
  TARGET_LINK_LIBRARIES( ${libname} Identifier)
ENDIF()

#------------------------------------------------------------------------------
# install
#------------------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

