#---------------------------------------------------------------------
# BeanUser is the place for user analysis + examples
#---------------------------------------------------------------------

# explicit listing of files
SET( CXX_FILES
  # List of test functions
  User1.cxx
  UserTest.cxx
  TestDb.cxx
  TestPID.cxx
  MagField.cxx
  JobInfo.cxx
  Rhopi.cxx
  EntryList.cxx
  TestEventTag.cxx
  TestAbsCor.cxx
  Etaetagamma.cxx
  #    the analysis of Jpsi -> rho pi  written for BEAN;
  #    analogue of the Boss: Analysis/Physics/RhopiAlg
  Rhopi.cxx
  #
  # -- local and old files --
  #     Semileptonic.cxx
  #     DoubleDTag.cxx
  #     Bhabha.cxx
  #     RadBhabha.cxx
  #
  #     SelectOmegaEta.cxx
  #     SelectOmegapi0.cxx
  #     Eta_TwoPi.cxx
  #     PipPimPi0Eta.cxx
  #     SelectJpsiX.cxx
  #
  SelectKKgg.cxx
  PsipJpsiPhiEta.cxx
  PsipPiPiKK.cxx
  PsipJpsiGammaEta.cxx
  PsipHelixParCor.cxx
  DelClonedTrks.cxx
  )

#---------------------------------------------------------------------
SET( libname BeanUser )

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${ROOTEVENTDATA_DIR}
  ${BeanCore_DIR}
  ${ParticleID_DIR}
  ${MagneticField_DIR}
  ${VertexFit_DIR}
  ${AbsCor_DIR}
  ${DatabaseSvc_DIR}
  ${EventTag_DIR}
  ${RscanDQ_DIR}
  ${TrackCorrection_DIR}
  ${ROOT_INCLUDE_DIR}
  ${CLHEP_INCLUDE_DIR}
  )

ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )

TARGET_LINK_LIBRARIES( ${libname}
  RootEventData
  BeanCore
  ParticleID
  MagneticField
  VertexFit
  AbsCor
  DatabaseSvc
  EventTag
  RscanDQ
  TrackCorrection
  ${ROOT_LIBRARIES}
  ${CLHEP_LIBRARIES}
  )

MAC_LIB_RPATH( ${libname} )

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

