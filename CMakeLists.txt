#------------------------------------------------------------------------------
# CMake configuration file to build Bean:
# see README.md how to use it
#------------------------------------------------------------------------------

CMAKE_MINIMUM_REQUIRED( VERSION 2.8.12 FATAL_ERROR )

# CMAKE_HOST_SYSTEM_NAME is about the system where cmake is being run
# CMAKE_SYSTEM_NAME is about the platform the script is building for
MESSAGE( STATUS
  "+++ BEAN: You are working on ${CMAKE_HOST_SYSTEM_NAME} system\n   "
  "+++ The program was tested under Linux, Windows and Darwin\n"
)

IF( CMAKE_HOST_SYSTEM_NAME MATCHES Darwin ) # MacOS
  # see "cmake/MacRpathPatch.cmake"
  SET ( CMAKE_MACOSX_RPATH ON )

  # this is for old versions of cmake, it is depricated now
  # MACOSX_RPATH (see cmake --help-policy CMP0042)
  #  CMAKE_POLICY( SET CMP0042 OLD )
  # LOCATION property for build targets (cmake --help-policy CMP0026)
  #  CMAKE_POLICY( SET CMP0026 OLD )
ENDIF()

PROJECT( BEAN )

#------------------------------------------------------------------------------
# Configuration parameters: use ccmake to change them
#------------------------------------------------------------------------------

# installation directory:
IF( NOT CMAKE_BUILD_TYPE )
  SET( CMAKE_INSTALL_PREFIX
     ${CMAKE_SOURCE_DIR}/workdir
     CACHE PATH
     "Choose the directory to install program"
     FORCE
  )
ENDIF()

# the full list of possible values of this flag are
# "DebugFull", "Debug", "Profile", "RelWithDebInfo" and "Release".
IF( NOT CMAKE_BUILD_TYPE )
   SET( CMAKE_BUILD_TYPE
      RelWithDebInfo
      CACHE STRING
      "Choose the type of build: RelWithDebInfo Release Debug"
      FORCE
   )
ENDIF()

# BOSS_VERSION is the version of BOSS program which was used
# to produce DST-files; the string BOSS_VERSION should contain
# three digits separated by dots; please do not specify
# the patch version of BOSS production.
# It defines:
# 1) the ROOTEVENTDATA_VERSION
# 2) versions of the algorithms in the Analysis folder;
#    (by setting the preprocessor variable BOSS_VER)
# Tested values: 6.6.3; 6.6.4; 6.6.5; 7.0.3; 7.0.4
SET( BOSS_VERSION
   6.6.4
   CACHE STRING
   "Version of BOSS in dst production"
)

OPTION( BUILTIN_SQLITE "Use included sqlite or system sqlite library" ON )
IF( NOT CMAKE_SYSTEM_NAME MATCHES Windows ) # not WINDOWS
  OPTION( ENABLE_CLHEP_XLINK "Enable RPATH-based link for CLHEP library" ON )
  OPTION( ENABLE_ROOT_XLINK "Enable RPATH-based link for ROOT libraries" ON )
ENDIF()

#------------------------------------------------------------------------------
# Set ROOTEVENTDATA_VERSION and BOSS_VER
#------------------------------------------------------------------------------
STRING( REGEX REPLACE "^([0-9]+)[.][0-9]+[.][0-9]+.*$" "\\1"
  BOSS_VER_MAJOR "${BOSS_VERSION}"
)
STRING( REGEX REPLACE "^[0-9][.]([0-9]+)[.][0-9]+.*$" "\\1"
  BOSS_VER_MINOR "${BOSS_VERSION}"
)
STRING( REGEX REPLACE "^[0-9]+[.][0-9]+[.]([0-9]+).*$" "\\1"
  BOSS_VER_MICRO "${BOSS_VERSION}"
)
SET( BOSS_VER "${BOSS_VER_MAJOR}${BOSS_VER_MINOR}${BOSS_VER_MICRO}" )
MESSAGE( STATUS "+++ BOSS_VER is ${BOSS_VER}" )

# ROOTEVENTDATA_VERSION defines format of dst-file:
#   BOSS_VERSION                        ROOTEVENTDATA_VERSION
#   6.5.2; 6.5.3; 6.5.5; 6.6.1          = BOSS_VERSION
#   6.6.2                               6.6.1
#   >= 6.6.3 < 7.0.5                    6.6.3
#   7.0.5                               7.0.5
#
SET( ROOTEVENTDATA_VERSION "6.6.3" )
IF( BOSS_VERSION  VERSION_EQUAL  "7.0.5" OR
    BOSS_VERSION  VERSION_LESS  "6.6.2" )
  SET( ROOTEVENTDATA_VERSION ${BOSS_VERSION} )
ELSEIF( BOSS_VERSION  VERSION_EQUAL  "6.6.2" )
  SET( ROOTEVENTDATA_VERSION "6.6.1" )
ENDIF()
#  MESSAGE( STATUS "+++ ROOTEVENTDATA_VERSION is ${ROOTEVENTDATA_VERSION}" )

SET( ROOTEVENTDATA RootEventData_${ROOTEVENTDATA_VERSION} )
IF( NOT EXISTS "${CMAKE_SOURCE_DIR}/${ROOTEVENTDATA}" )
   MESSAGE( FATAL_ERROR
   "\"${ROOTEVENTDATA}\" does not exist in folder:\n"
   "\"${CMAKE_SOURCE_DIR}\"\n"
   "The BOSS_VERSION=${BOSS_VERSION} is not supported or "
   "is specified incorrectly.\n" )
ENDIF()
MESSAGE( STATUS "+++ ROOTEVENTDATA: ../${ROOTEVENTDATA}" )


#------------------------------------------------------------------------------
# Try to find ROOT and CLHEP
#------------------------------------------------------------------------------
SET( CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake )

SET( ROOT_CONFIG_SEARCHPATH
   ""
   CACHE PATH
   "Path to directory with root-config program"
)
# uncomment for detailed printing of found root-configuration
# SET( FIND_ROOT_VERBOSE "VERBOSE" )
INCLUDE( FindROOT )

SET( CLHEP_SEARCHPATH
   "$ENV{CLHEP_DIR}"
   CACHE PATH
   "Path to directory where clhep is installed"
)
INCLUDE( FindCLHEP )

#------------------------------------------------------------------------------
# Compilation and link/libraries options
#------------------------------------------------------------------------------

# A minimum CXX version:

IF( CMAKE_COMPILER_IS_GNUCC )
  # Since CMake 2.8.10 there is the CMAKE_CXX_COMPILER_VERSION variable
  # but keep it for cmake 2.6
  EXECUTE_PROCESS( COMMAND ${CMAKE_C_COMPILER} -dumpversion
                   OUTPUT_VARIABLE GCC_VERSION )
#  MESSAGE( STATUS "gcc_verion= ${GCC_VERSION}" )
# c++14 mode is the default in GCC 6.1 up until GCC 10 (including)
IF( GCC_VERSION VERSION_LESS 6.1 )
  IF( GCC_VERSION VERSION_GREATER 4.7 )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
  ELSEIF( GCC_VERSION VERSION_GREATER 4.3 )
    SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x" )
  ELSEIF()
    MESSAGE( FATAL_ERROR "Insufficient gcc version: ${GCC_VERSION}"
             " we use C++11 features" )
  ENDIF()
ENDIF()

  # just to see all compilation warnings:
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
ENDIF()

# for Clang or AppleClang
IF ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" )
  SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++11" )
ENDIF()

# WINDOWS specific
IF( CMAKE_SYSTEM_NAME MATCHES Windows )
  # _USE_MATH_DEFINES - to use Math Constants
  # _CRT_SECURE_NO_WARNINGS - to eliminate deprecation warnings for the older,
  #                           less secure functions
  ADD_DEFINITIONS( -D_USE_MATH_DEFINES -D_CRT_SECURE_NO_WARNINGS )
  IF( CMAKE_BUILD_TYPE MATCHES Release )
    ADD_DEFINITIONS( -MD )
  ELSE()
    ADD_DEFINITIONS( -MDd )
  ENDIF()
ENDIF()

# Set the common place where all libraries to be built.
# Please do not to be confused it with the place where
# the libraries will be installed later.
SET( LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib )

# Set actual place where the libraries will be installed
SET( install_lib_dir "BeanLib_${BOSS_VERSION}" )

IF( NOT CMAKE_SYSTEM_NAME MATCHES Windows ) # not WINDOWS
  SET( LD_LIBRARY_PATH "LD_LIBRARY_PATH" )
  IF( CMAKE_SYSTEM_NAME MATCHES Darwin ) # MacOS
    SET( LD_LIBRARY_PATH "DYLD_LIBRARY_PATH" )
  ENDIF()

  # generate setup.sh
  CONFIGURE_FILE(
     ${CMAKE_SOURCE_DIR}/cmake/setup.sh.cmake
     ${PROJECT_BINARY_DIR}/setup_${BOSS_VERSION}.sh
     @ONLY
  )
ENDIF()

#------------------------------------------------------------------------------
# RPATH settings
#------------------------------------------------------------------------------
INCLUDE( MacRpathPatch )

IF( NOT CMAKE_SYSTEM_NAME MATCHES Windows ) # not WINDOWS
  IF( ENABLE_ROOT_XLINK )
     SET( ROOTXLINK "-Xlinker -rpath -Xlinker ${ROOT_LIBRARY_DIR}" )
     SET( CMAKE_SHARED_LINKER_FLAGS
          "${CMAKE_SHARED_LINKER_FLAGS} ${ROOTXLINK}" )
  ENDIF()
  IF( ENABLE_CLHEP_XLINK )
     SET( CLHEPXLINK "-Xlinker -rpath -Xlinker ${CLHEP_LIBRARY_DIR}" )
     SET( CMAKE_SHARED_LINKER_FLAGS
          "${CMAKE_SHARED_LINKER_FLAGS} ${CLHEPXLINK}" )
  ENDIF()

  SET( CMAKE_SKIP_BUILD_RPATH FALSE )

  # append directories in the linker search path
  SET( CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE )

  SET( CMAKE_INSTALL_RPATH
       "${CMAKE_INSTALL_PREFIX}/${install_lib_dir}" )
ENDIF()

#------------------------------------------------------------------------------
# build libraries
#------------------------------------------------------------------------------
ADD_DEFINITIONS( -DBEAN )
ADD_DEFINITIONS( -DBOSS_VER=${BOSS_VER} )

ADD_SUBDIRECTORY( ${ROOTEVENTDATA} )
SET( ROOTEVENTDATA_DIR  ${CMAKE_SOURCE_DIR}/${ROOTEVENTDATA} )

ADD_SUBDIRECTORY( BeanCore )
SET( BeanCore_DIR  ${CMAKE_SOURCE_DIR}/BeanCore )

# -- Analysis --
IF( BOSS_VERSION  VERSION_GREATER  "7.0.0" )
  ADD_SUBDIRECTORY( Analysis/Identifier )
  SET( Identifier_DIR  ${CMAKE_SOURCE_DIR}/Analysis/Identifier )
ENDIF()

ADD_SUBDIRECTORY( Analysis/AbsCor )
SET( AbsCor_DIR  ${CMAKE_SOURCE_DIR}/Analysis/AbsCor )

ADD_SUBDIRECTORY( Analysis/DatabaseSvc )
SET( DatabaseSvc_DIR  ${CMAKE_SOURCE_DIR}/Analysis/DatabaseSvc )

ADD_SUBDIRECTORY( Analysis/MagneticField )
SET( MagneticField_DIR  ${CMAKE_SOURCE_DIR}/Analysis/MagneticField )

ADD_SUBDIRECTORY( Analysis/EventTag )
SET( EventTag_DIR  ${CMAKE_SOURCE_DIR}/Analysis/EventTag )

IF( BOSS_VERSION  VERSION_GREATER  "7.0.0" )
  ADD_SUBDIRECTORY( Analysis/ParticleID )
  SET( ParticleID_DIR  ${CMAKE_SOURCE_DIR}/Analysis/ParticleID )
ELSE()
  ADD_SUBDIRECTORY( Analysis/ParticleID_boss6 )
  SET( ParticleID_DIR  ${CMAKE_SOURCE_DIR}/Analysis/ParticleID_boss6 )
ENDIF()

ADD_SUBDIRECTORY( Analysis/VertexFit )
SET( VertexFit_DIR  ${CMAKE_SOURCE_DIR}/Analysis/VertexFit )

ADD_SUBDIRECTORY( Analysis/RscanDQ )
SET( RscanDQ_DIR  ${CMAKE_SOURCE_DIR}/Analysis/RscanDQ )

ADD_SUBDIRECTORY( Analysis/TrackCorrection )
SET( TrackCorrection_DIR  ${CMAKE_SOURCE_DIR}/Analysis/TrackCorrection )

# -- User --
ADD_SUBDIRECTORY( BeanUser )

#------------------------------------------------------------------------------
# BEAN executable
#------------------------------------------------------------------------------

# program name
SET( prog bean_${BOSS_VERSION}.exe )
IF( CMAKE_SYSTEM_NAME MATCHES Windows )
  # in windows the `exe` extention is added by default
  SET( prog bean_${BOSS_VERSION} )
ENDIF()

INCLUDE_DIRECTORIES(
   ${ROOTEVENTDATA_DIR}
   ${BeanCore_DIR}
   ${ROOT_INCLUDE_DIR}
   ${CLHEP_INCLUDE_DIR}
   ${DatabaseSvc_DIR}
)

SET( MAIN_SRCS "BeanCore/main.cxx" )

# add getopt() port for MSVC:
# http://www.codeproject.com/KB/cpp/getopt4win.aspx
IF( CMAKE_SYSTEM_NAME MATCHES Windows )
  ADD_DEFINITIONS( -DSTATIC_GETOPT )
  LIST( APPEND MAIN_SRCS "BeanCore/win_getopt/getopt.c" )
ENDIF()

ADD_DEFINITIONS( -DBEANBASE="${CMAKE_SOURCE_DIR}" )

ADD_EXECUTABLE( ${prog} ${MAIN_SRCS} )

IF( CMAKE_SYSTEM_NAME MATCHES Linux )
  # is it platform or compiler specific?
  SET( CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-as-needed" )
ENDIF()

TARGET_LINK_LIBRARIES( ${prog}
   RootEventData
   BeanCore
   BeanUser
   AbsCor
   DatabaseSvc
   MagneticField
   EventTag
   ParticleID
   VertexFit
   ${ROOT_LIBRARIES}
   ${CLHEP_LIBRARIES}
)

IF( NOT CMAKE_SYSTEM_NAME MATCHES Windows ) # not WINDOWS
  SET_TARGET_PROPERTIES( ${prog}
     PROPERTIES  LINK_FLAGS  "${ROOTXLINK} ${CLHEPXLINK}"
  )
  MAC_PROG_RPATH( ${prog} )
ENDIF()

# final message (must be simple!) after compilation is done
ADD_CUSTOM_COMMAND( TARGET ${prog} POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E cmake_echo_color --red --bold
  "Use   make install   to install the program!"
  COMMENT "Compilation is over!"
)

#------------------------------------------------------------------------------
# install
#------------------------------------------------------------------------------
INSTALL( DIRECTORY DESTINATION ${install_lib_dir} )
INSTALL( TARGETS ${prog}  DESTINATION ${CMAKE_INSTALL_PREFIX} )
INSTALL( FILES ${PROJECT_BINARY_DIR}/setup_${BOSS_VERSION}.sh
         DESTINATION ${CMAKE_INSTALL_PREFIX} )

ADD_CUSTOM_TARGET( install_bean
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake"
)

#------------------------------------------------------------------------------
# PROOF
# Note: the name of par-file MUST match the name of the directory
#------------------------------------------------------------------------------

# PROOF_A_FILES is the list of data files that could be used
# in the BEAN from the directory "Analysis"

FILE( GLOB PROOF_A_FILES
   RELATIVE "${CMAKE_SOURCE_DIR}" "Analysis/PROOF-INF/*.*" )
FILE( GLOB TMP_FILES
   RELATIVE "${CMAKE_SOURCE_DIR}" "Analysis/DatabaseSvc/dat/*.db" )
LIST( APPEND PROOF_A_FILES ${TMP_FILES} )
FILE( GLOB TMP_FILES
   RELATIVE "${CMAKE_SOURCE_DIR}" "Analysis/MagneticField/dat/*.dat" )
LIST( APPEND PROOF_A_FILES ${TMP_FILES} )
FILE( GLOB TMP_FILES
   RELATIVE "${CMAKE_SOURCE_DIR}" "Analysis/MagneticField/dat/*.txt" )
LIST( APPEND PROOF_A_FILES ${TMP_FILES} )
FILE( GLOB_RECURSE TMP_FILES
   RELATIVE "${CMAKE_SOURCE_DIR}" "Analysis/ParticleID/share/*.txt" )
LIST( APPEND PROOF_A_FILES ${TMP_FILES} )
FILE( GLOB TMP_FILES
   RELATIVE "${CMAKE_SOURCE_DIR}" "Analysis/AbsCor/dat/*.txt" )
LIST( APPEND PROOF_A_FILES ${TMP_FILES} )
FILE( GLOB TMP_FILES
   RELATIVE "${CMAKE_SOURCE_DIR}" "Analysis/EventTag/share/*.codes" )
LIST( APPEND PROOF_A_FILES ${TMP_FILES} )
LIST( APPEND PROOF_A_FILES "Analysis/EventTag/share/pdt.table" )

# MESSAGE( STATUS "PROOF_A_FILES = ${PROOF_A_FILES}" )

SET( ANALYSIS_PAR "${CMAKE_INSTALL_PREFIX}/par/Analysis.par" )
ADD_CUSTOM_COMMAND(
   OUTPUT   ${ANALYSIS_PAR}
   COMMAND "${CMAKE_COMMAND}" -E make_directory ${CMAKE_INSTALL_PREFIX}/par
   COMMAND "${CMAKE_COMMAND}" -E tar -czf ${ANALYSIS_PAR} ${PROOF_A_FILES}
   DEPENDS  ${PROOF_A_FILES}
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
   COMMENT "Creating ${ANALYSIS_PAR}"
)

# PROOF_B_FILES is the list of data files that could be used
# in the BEAN from the directory "BeanUser"

FILE( GLOB PROOF_B_FILES
   RELATIVE "${CMAKE_SOURCE_DIR}" "BeanUser/PROOF-INF/*.*" )
LIST( APPEND PROOF_B_FILES "BeanUser/mydecay.codes" )

SET( BEANUSER_PAR "${CMAKE_INSTALL_PREFIX}/par/BeanUser.par" )
ADD_CUSTOM_COMMAND(
   OUTPUT   ${BEANUSER_PAR}
   COMMAND "${CMAKE_COMMAND}" -E make_directory ${CMAKE_INSTALL_PREFIX}/par
   COMMAND "${CMAKE_COMMAND}" -E tar -czf ${BEANUSER_PAR} ${PROOF_B_FILES}
   DEPENDS  ${PROOF_B_FILES}
   WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
   COMMENT "Creating ${BEANUSER_PAR}"
)

SET( LIB_PAR "${CMAKE_INSTALL_PREFIX}/par/${install_lib_dir}.par" )
SET( PROOF_INF_LIB "${CMAKE_INSTALL_PREFIX}/${install_lib_dir}/PROOF-INF" )
ADD_CUSTOM_COMMAND(
   OUTPUT   ${PROOF_INF_LIB}/BUILD.sh ${PROOF_INF_LIB}/SETUP.C
   COMMAND "${CMAKE_COMMAND}" -E make_directory ${PROOF_INF_LIB}
   COMMAND "${CMAKE_COMMAND}" -E
           copy BUILD_lib.sh ${PROOF_INF_LIB}/BUILD.sh
   COMMAND "${CMAKE_COMMAND}" -E
           copy SETUP_lib.C ${PROOF_INF_LIB}/SETUP.C
   DEPENDS ${CMAKE_SOURCE_DIR}/cmake/BUILD_lib.sh
           ${CMAKE_SOURCE_DIR}/cmake/SETUP_lib.C
   WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}/cmake"
   COMMENT "Creating ${PROOF_INF_LIB}"
)

ADD_CUSTOM_COMMAND(
   OUTPUT   ${LIB_PAR}
   COMMAND "${CMAKE_COMMAND}" -E make_directory ${CMAKE_INSTALL_PREFIX}/par
   COMMAND "${CMAKE_COMMAND}" -E tar -czf ${LIB_PAR} ${install_lib_dir}
   DEPENDS ${PROOF_INF_LIB}/BUILD.sh ${PROOF_INF_LIB}/SETUP.C
   WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}
   COMMENT "Creating ${LIB_PAR}"
)

ADD_CUSTOM_TARGET( proofbean
   DEPENDS ${ANALYSIS_PAR} ${BEANUSER_PAR} ${LIB_PAR}
)
ADD_DEPENDENCIES( proofbean install_bean )

#------------------------------------------------------------------------------
# uninstall -> delete "by hand"
#------------------------------------------------------------------------------
ADD_CUSTOM_TARGET( uninstall
   COMMAND "${CMAKE_COMMAND}" -E remove -f ${prog}
   COMMAND "${CMAKE_COMMAND}" -E remove -f setup_${BOSS_VERSION}.sh
   COMMAND "${CMAKE_COMMAND}" -E remove_directory ${install_lib_dir}
   COMMAND "${CMAKE_COMMAND}" -E remove_directory par
   WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}
)
