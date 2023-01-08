#---------------------------------------------------------------------
# RscanDQ:
#  https://docbes3.ihep.ac.cn/~tauqcdgroup/index.php/Data_Samples
#  When you using the 2015 Rscan data, you must use this package
#  to exclude bad runs :
#  "/besfs/users/liud/workarea/workarea665p01/RscanDQ/RscanDQ-00-00-01"
#---------------------------------------------------------------------

SET( libname RscanDQ )

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  )

#---------------------------------------------------------------------
# the list of source files
FILE( GLOB CXX_FILES src/*.cxx )

#---------------------------------------------------------------------
#
IF( CMAKE_SYSTEM_NAME MATCHES Windows ) # WINDOWS
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
ELSE()
  ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )
  MAC_LIB_RPATH( ${libname} )
ENDIF()

#---------------------------------------------------------------------
# install
#---------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

