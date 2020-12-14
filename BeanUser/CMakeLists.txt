#--------------------------------------------------------------------------
# BeanUser is the place for user analysis + examples
#--------------------------------------------------------------------------

# explicit listing of files
SET( CXX_FILES
   User1.cxx
   UserTest.cxx
   TestDb.cxx
   TestPID.cxx
   MagField.cxx
   JobInfo.cxx
   Rhopi.cxx
   Bhabha.cxx
   RadBhabha.cxx
   EntryList
   TestEventTag.cxx
   Etaetagamma.cxx
#
#   Semileptonic.cxx
#   DoubleDTag.cxx
#
#   see ~/study/bes3/ARCHIVE/BeanUser/
#   Select3pi.cxx
#   SelectK3pi.cxx
#   SelectK3pi_2VF.cxx
#   SelectKsKpi.cxx
#   SelectPhiEta958.cxx
#   SelectOmegaEta958.cxx
#   SelectKKpi0.cxx
#
#     SelectOmegapi0.cxx
#     SelectOmegaEta.cxx
   SelectKKgg.cxx

   PsipJpsiPhiEta.cxx
   PsipPiPiKK.cxx
   PsipJpsiGammaEta.cxx
#     PsipMuMuGamma.cxx
#
   Eta_TwoPi.cxx
   PipPimPi0Eta.cxx
#     SelectJpsiX.cxx
)

#--------------------------------------------------------------------------
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
   ${ROOT_LIBRARIES}
   ${CLHEP_LIBRARIES}
)

MAC_LIB_RPATH( ${libname} )

#--------------------------------------------------------------------------
# install
#--------------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

