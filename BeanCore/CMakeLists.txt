#---------------------------------------------------------------------
# BeanCore: library provides the basic functionality of the bean
#           and contains transitional classes Root->Analysis
#---------------------------------------------------------------------

SET( libname BeanCore )

SET( INCLUDE_DIRECTORIES
   ${CMAKE_CURRENT_SOURCE_DIR}
   ${ROOTEVENTDATA_DIR}
   ${ROOT_INCLUDE_DIR}
   ${CLHEP_INCLUDE_DIR}
)
INCLUDE_DIRECTORIES( ${INCLUDE_DIRECTORIES} )

# MESSAGE( " +++Bean INCLUDE_DIRECTORIES= ${INCLUDE_DIRECTORIES}" )

#---------------------------------------------------------------------
# the list of headers are needed to generate the root dictionary
SET( LINKDEFHDRS        DstFormat.h ReadDst.h Bean.h
                        MasterInfoHack.h TMergeableMap.h )

IF( ROOT_MAJOR_VERS LESS 6 )
  # the name of the dictionary
  SET( DICTIONARY ${libname}_rootcint.cxx )
  ROOT_GENERATE_DICTIONARY(
     ${DICTIONARY}
     "${LINKDEFHDRS}"
     "${INCLUDE_DIRECTORIES}"
     ""
  )
ELSE()
  ROOT6_GENERATE_DICTIONARY(
     ${libname}
     "${LINKDEFHDRS}"
     "${INCLUDE_DIRECTORIES}"
     ""
  )
  SET( DICTIONARY ${libname}_rdict.cxx )
  SET( DictNamePcm ${libname}_rdict.pcm )
ENDIF()

#---------------------------------------------------------------------
# the list of source files
SET( CXX_FILES  DstFormat.cxx ReadDst.cxx Bean.cxx
                MasterInfoHack.cxx TMergeableMap.cxx
                DstEvtRecTracks.cxx TofHitStatus.cxx
                BMdcTrack.cxx BMdcKalTrack.cxx BEmcTrack.cxx
                BExtTrack.cxx
                unix_cp.cxx )

# add the dictionary to the list of source files
LIST( APPEND CXX_FILES ${DICTIONARY} )

#---------------------------------------------------------------------
#
IF( CMAKE_SYSTEM_NAME MATCHES Windows ) # WINDOWS
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
  TARGET_LINK_LIBRARIES( ${libname}
    RootEventData
    ${ROOT_LIBRARIES}
    ${CLHEP_LIBRARIES}
  )
ELSE()
  ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )
  SET( Bean_ROOT_LIBRARIES "${ROOT_LIBRARIES} -lProofPlayer -lProof" )
  TARGET_LINK_LIBRARIES( ${libname}
    RootEventData
    ${Bean_ROOT_LIBRARIES}
    ${CLHEP_LIBRARIES}
  )
  MAC_LIB_RPATH( ${libname} )
ENDIF()

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

IF( ROOT_MAJOR_VERS GREATER 5 )
  INSTALL( FILES ${CMAKE_CURRENT_BINARY_DIR}/${DictNamePcm}
           DESTINATION ${install_lib_dir} )
ENDIF()
