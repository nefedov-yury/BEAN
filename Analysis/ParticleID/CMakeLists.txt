#---------------------------------------------------------------------
# ParticleID: library for particle identification
#              <=> boss/Analysis/ParticleID
#---------------------------------------------------------------------

SET( libname ParticleID )

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOTEVENTDATA_DIR}
  ${BeanCore_DIR}
  ${CLHEP_INCLUDE_DIR}
  ${ROOT_INCLUDE_DIR}
  )

#---------------------------------------------------------------------
# the list of source files
FILE( GLOB CXX_FILES src/*.cxx )

#---------------------------------------------------------------------
#
IF( CMAKE_SYSTEM_NAME MATCHES Windows ) # WINDOWS
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
  TARGET_LINK_LIBRARIES( ${libname}
    RootEventData
    BeanCore
    ${ROOT_LIBRARIES}
    ${CLHEP_LIBRARIES}
    )
ELSE()
  ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )
  SET( ParticleID_ROOT_LIBRARIES "${ROOT_LIBRARIES} -lMLP" )
  TARGET_LINK_LIBRARIES( ${libname}
    RootEventData
    BeanCore
    ${ParticleID_ROOT_LIBRARIES}
    ${CLHEP_LIBRARIES}
    )
  MAC_LIB_RPATH( ${libname} )
ENDIF()

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

