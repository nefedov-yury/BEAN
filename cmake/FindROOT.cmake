###########################################################################
# - Try to find ROOT installation
#
# For UNIX-like system this module tries to find the root-config script.
# Then, using this script extracts the meta-information about ROOT
# installation.
#
# In Windows we rely on environment variable ROOTSYS and
# set the directories bin, lib and include with respect to
# the folder ${ROOTSYS), if it exists.
# (The root-config script requires bash and a number of standard
#  UNIX tools like test, sed, awk and so on)
#
# Input:
#       The variable ROOT_CONFIG_SEARCHPATH may contain the full path to
#       root-config script. Otherwise, we check directory $ROOTSYS/bin
#       using the environment variable ROOTSYS.
#       FIND_ROOT_VERBOSE - set verbose output for this script
#
# Output:
#       ROOTSYS  - system has the ROOT and $ROOTSYS is
#                  the path to installation directory
#
#       ROOT_CONFIG_EXECUTABLE - root-config with full path to it
#       ROOT_CINT_EXECUTABLE   - rootcint (or rootcling for ROOT-6)
#
#       ROOT_LIBRARY_DIR - the library directory
#       ROOT_BINARY_DIR  - the executable directory
#       ROOT_INCLUDE_DIR - the header directory
#       ROOT_LIBRARIES   - regular ROOT libraries
#       ROOT_GLIBS       - regular + GUI ROOT libraries
#       ROOT_CFLAGS      - extra compiler flags
#
#       ROOTVERSION - the version of ROOT (something like 5.14/00h) and
#       ROOT_MAJOR_VERS,ROOT_MINOR_VERS,ROOT_PATCH_VERS
#                   - the corresponding components of the version string
#
#       ROOT_GENERATE_DICTIONARY -  function for building ROOT dictionary
#
# Nefedov dot Yury at jinr dot ru
# Based on a script written by the F.Uhlig@gsi.de (fairroot.gsi.de)
###########################################################################

# if we already found ROOT do nothing
IF( NOT ROOTSYS )
  MESSAGE( STATUS "Looking for Root..." )

  # function we will use
  FUNCTION( _root_check  File )
    IF( NOT EXISTS "${File}" )
      SET( ROOTSYS, NOTFOUND )
      MESSAGE( FATAL_ERROR "${File} does not exist. "
        "ROOT must not be installed correctly."
      )
    ENDIF()
  ENDFUNCTION( _root_check )


  IF( CMAKE_SYSTEM_NAME MATCHES Windows ) # WINDOWS
    IF( MSVC )
      MESSAGE( STATUS "Compiler: MSVC, version: " ${MSVC_VERSION} )
    ELSE()
      MESSAGE( FATAL_ERROR "Only MSVC is supported on Windows" )
    ENDIF()

    # Set ROOTSYS
    IF( ROOT_CONFIG_SEARCHPATH )
      STRING( REGEX REPLACE "(^.*)[/\]bin" "\\1"
        ROOTSYS ${ROOT_CONFIG_SEARCHPATH}
      )
    ELSE()
      SET( ROOTSYS $ENV{ROOTSYS} )
    ENDIF()
    FILE( TO_CMAKE_PATH "${ROOTSYS}" ROOTSYS )
    MESSAGE( STATUS "Looking for Root... ROOTSYS=\n\t${ROOTSYS}" )

    IF( ROOTSYS )
      _root_check( ${ROOTSYS} )
      # SET( $ENV{ROOTSYS} ${ROOTSYS} ) # set path ??
    ELSE()
      MESSAGE( FATAL_ERROR "The environment variable ROOTSYS"
        " _must be_ defined on Windows system"
      )
    ENDIF()

    # Set ROOT_BINARY_DIR
    SET( ROOT_BINARY_DIR "${ROOTSYS}/bin" )
    _root_check( ${ROOT_BINARY_DIR} )

    # Set ROOT_LIBRARY_DIR
    SET( ROOT_LIBRARY_DIR "${ROOTSYS}/lib" )
    _root_check( ${ROOT_LIBRARY_DIR} )

    # Set ROOT_INCLUDE_DIR
    SET( ROOT_INCLUDE_DIR "${ROOTSYS}/include" )
    _root_check( ${ROOT_INCLUDE_DIR} )

    # Check file RVersion.h and set ROOTVERSION
    _root_check( "${ROOT_INCLUDE_DIR}/RVersion.h" )
    FILE( READ "${ROOT_INCLUDE_DIR}/RVersion.h" contents )
    #  MESSAGE(STATUS "+++ contents=${contents} +++" )
    STRING(
      REGEX REPLACE ".*define[ ]+ROOT_RELEASE[ ]+[\"]([^\"]*).*" "\\1"
      ROOTVERSION  "${contents}"
    )

    # Set ROOT_LIBRARIES and ROOT_GLIBS == ROOT_LIBRARIES
    FILE( GLOB ROOT_LIBRARIES ${ROOT_LIBRARY_DIR}/*.lib )
    SET( ROOT_GLIBS ${ROOT_LIBRARIES} )

    # Set ROOT_CFLAGS
    SET( ROOT_CFLAGS "-I${ROOT_INCLUDE_DIR}" )

  ELSE() # LINUX and MACOS

    IF( NOT ROOT_CONFIG_SEARCHPATH )
      SET( ROOT_CONFIG_SEARCHPATH $ENV{ROOTSYS}/bin )
    ENDIF()

    FIND_PROGRAM( ROOT_CONFIG_EXECUTABLE
      NAMES root-config
      PATHS ${ROOT_CONFIG_SEARCHPATH}
      NO_DEFAULT_PATH
    )

    IF( ROOT_CONFIG_EXECUTABLE )
      MESSAGE( STATUS
        "Looking for Root... found root-config:\n\t"
        " ${ROOT_CONFIG_EXECUTABLE}"
      )

      # Set ROOTSYS
      STRING( REGEX REPLACE "(^.*)/bin/root-config" "\\1"
        ROOTSYS ${ROOT_CONFIG_EXECUTABLE}
      )

      FUNCTION( _root_config ARG OUT )
        #  EXEC_PROGRAM( ${ROOT_CONFIG_EXECUTABLE} ARGS ${ARG}
        EXECUTE_PROCESS( COMMAND ${ROOT_CONFIG_EXECUTABLE} ${ARG}
          OUTPUT_VARIABLE tmp
          OUTPUT_STRIP_TRAILING_WHITESPACE )
        IF( CMAKE_SYSTEM_NAME MATCHES Windows )
          STRING(REPLACE "\r\n" "" tmp "${tmp}")
        ELSE()
          STRING(REPLACE "\n" "" tmp "${tmp}")
        ENDIF()
        SET( ${OUT} ${tmp} PARENT_SCOPE )
        #  MESSAGE(STATUS "_root_config ${ARG} => ${OUT} = ${tmp}" )
      ENDFUNCTION( _root_config )

      # Set ROOTVERSION
      _root_config( --version ROOTVERSION )

      # Set ROOT_LIBRARY_DIR
      _root_config( --libdir ROOT_LIBRARY_DIR )
      _root_check( ${ROOT_LIBRARY_DIR} )

      # Set ROOT_BINARY_DIR
      _root_config( --bindir ROOT_BINARY_DIR )
      _root_check( ${ROOT_BINARY_DIR} )

      # Set ROOT_INCLUDE_DIR
      _root_config( --incdir ROOT_INCLUDE_DIR )
      _root_check( ${ROOT_INCLUDE_DIR} )

      # Set ROOT_LIBRARIES and ROOT_GLIBS
      _root_config( --libs ROOT_LIBRARIES )
      _root_config( --glibs ROOT_GLIBS )

      # Set ROOT_CFLAGS
      _root_config( --cflags ROOT_CFLAGS )
      #  MESSAGE( STATUS "ROOT_CFLAGS= ${ROOT_CFLAGS}" )

    ELSE ()
      SET( ROOTSYS, NOTFOUND )
      MESSAGE( FATAL_ERROR
        "ROOT not installed in the searchpath and ROOTSYS is not set. "
        "Please specify the ROOT_CONFIG_SEARCHPATH variable the path "
        "to the root-config program or set the environment variable ROOTSYS."
      )
    ENDIF( ROOT_CONFIG_EXECUTABLE )

  ENDIF() # OPERATION SYSTEM

  MESSAGE( STATUS
    "Looking for Root... found version: ${ROOTVERSION}" )

  # now parse the ROOTVERSION string into variables
  # ( note that two backslashes (\\) are required to get a backslash
  #   through argument parsing )
  STRING( REGEX REPLACE "^([0-9]+)[.][0-9]+[/][0-9]+.*$" "\\1"
    ROOT_MAJOR_VERS     "${ROOTVERSION}"
  )
  STRING( REGEX REPLACE "^[0-9]+[.]([0-9]+)[/][0-9]+.*$" "\\1"
    ROOT_MINOR_VERS     "${ROOTVERSION}"
  )
  STRING( REGEX REPLACE "^[0-9]+[.][0-9]+[/]([0-9]+).*$" "\\1"
    ROOT_PATCH_VERS     "${ROOTVERSION}"
  )

  # check the executables of rootcint and set ROOT_CINT_EXECUTABLE
  SET( _rootcint rootcint )
  IF( ROOT_MAJOR_VERS GREATER 5 )
    SET( _rootcint rootcling )
  ENDIF()
  FIND_PROGRAM( ROOT_CINT_EXECUTABLE
    NAMES ${_rootcint}
    PATHS ${ROOT_BINARY_DIR}
    NO_DEFAULT_PATH
  )
  MESSAGE( STATUS "Looking for Root... found ${_rootcint}:\n\t "
    "${ROOT_CINT_EXECUTABLE}"
  )

  IF( FIND_ROOT_VERBOSE )
    MESSAGE( STATUS "Looking for Root... \n"
      "\t version components: ${ROOT_MAJOR_VERS} ${ROOT_MINOR_VERS} "
      "${ROOT_PATCH_VERS}\n"
      "\t ROOT_LIBRARY_DIR= ${ROOT_LIBRARY_DIR}\n"
      "\t ROOT_BINARY_DIR= ${ROOT_BINARY_DIR}\n"
      "\t ROOT_INCLUDE_DIR= ${ROOT_INCLUDE_DIR}\n"
      "\t ROOT_LIBRARIES= ${ROOT_LIBRARIES}\n"
      "\t ROOT_GLIBS= ${ROOT_GLIBS}\n"
    )
  ENDIF()

ENDIF( NOT ROOTSYS )


###########################################################################
#
# function for building ROOT dictionaries
#  Parameters:
#       DictNameCxx  - the C++ name of CINT dictionary file
#       IncludeFiles - the list of include files
#       IncludeDirs  - the include file directories to be searched
#       LinkDefH     - the LinkDef.h file (optional)
#
#  Output files:  two dictionary files ${DictNameCxx} and
#                 complimentary header file with ".h" suffix
#
###########################################################################

FUNCTION( ROOT_GENERATE_DICTIONARY
  DictNameCxx IncludeFiles IncludeDirs LinkDefH )

  # add "-I" before each include file directory
  SET( INCLUDE_DIRS )
  FOREACH( dir ${IncludeDirs} )
    SET( INCLUDE_DIRS ${INCLUDE_DIRS} -I${dir} )
  ENDFOREACH()

  STRING( REGEX REPLACE "^(.*)[.](.*)$" "\\1.h"
    DictNameH "${DictNameCxx}" )
  SET( OUTFILES ${DictNameCxx} ${DictNameH} )

  IF( CMAKE_SYSTEM_NAME MATCHES Linux )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${OUTFILES}
      COMMAND LD_LIBRARY_PATH=${ROOT_LIBRARY_DIR} ROOTSYS=${ROOTSYS}
              ${ROOT_CINT_EXECUTABLE}
      ARGS -f ${DictNameCxx}
           -c ${INCLUDE_DIRS} ${IncludeFiles} ${LinkDefH}
      DEPENDS ${IncludeFiles} ${LinkDefH}
    )
  ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${OUTFILES}
      COMMAND ${ROOT_CINT_EXECUTABLE}
      ARGS -f ${DictNameCxx}
           -c ${INCLUDE_DIRS} ${IncludeFiles} ${LinkDefH}
      DEPENDS ${IncludeFiles} ${LinkDefH}
    )
  ELSEIF( CMAKE_SYSTEM_NAME MATCHES Darwin )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${OUTFILES}
      COMMAND DYLD_LIBRARY_PATH=${ROOT_LIBRARY_DIR} ROOTSYS=${ROOTSYS}
              ${ROOT_CINT_EXECUTABLE}
      ARGS -f ${DictNameCxx}
           -c ${INCLUDE_DIRS} ${IncludeFiles} ${LinkDefH}
      DEPENDS ${IncludeFiles} ${LinkDefH}
    )
  ENDIF()

ENDFUNCTION( ROOT_GENERATE_DICTIONARY )


###########################################################################
#
# function for building ROOT6 dictionaries
#  Parameters:
#       TagLibName   - the target library name
#       IncludeFiles - the list of include files
#       IncludeDirs  - the include file directories to be searched
#       LinkDefH     - the LinkDef.h file
#
#  Output files:  two dictionary files with names build from TagLibName:
#                 TagLibName_rdict.cxx & TagLibName_rdict.pcm
#
###########################################################################

FUNCTION( ROOT6_GENERATE_DICTIONARY
    TagLibName IncludeFiles IncludeDirs LinkDefH )

  # add "-I" before each include file directory
  SET( INCLUDE_DIRS )
  FOREACH( dir ${IncludeDirs} )
    SET( INCLUDE_DIRS ${INCLUDE_DIRS} -I${dir} )
  ENDFOREACH()

  SET( DictNameCxx ${TagLibName}_rdict.cxx )
  SET( DictNamePcm ${TagLibName}_rdict.pcm )
  SET( OUTFILES ${DictNameCxx} ${DictNamePcm} )

  IF( CMAKE_SYSTEM_NAME MATCHES Linux )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${OUTFILES}
      COMMAND LD_LIBRARY_PATH=${ROOT_LIBRARY_DIR} ROOTSYS=${ROOTSYS}
              ${ROOT_CINT_EXECUTABLE}
      ARGS -f ${DictNameCxx} -s ${TagLibName}
        ${INCLUDE_DIRS} ${IncludeFiles} ${LinkDefH}
      DEPENDS ${IncludeFiles} ${LinkDefH}
    )
  ELSEIF( CMAKE_SYSTEM_NAME MATCHES Windows )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${OUTFILES}
      COMMAND ${ROOT_CINT_EXECUTABLE}
      ARGS -f ${DictNameCxx} -s ${TagLibName}
        ${INCLUDE_DIRS} ${IncludeFiles} ${LinkDefH}
      DEPENDS ${IncludeFiles} ${LinkDefH}
    )
  ELSEIF( CMAKE_SYSTEM_NAME MATCHES Darwin )
    ADD_CUSTOM_COMMAND(
      OUTPUT ${OUTFILES}
      COMMAND DYLD_LIBRARY_PATH=${ROOT_LIBRARY_DIR} ROOTSYS=${ROOTSYS}
              ${ROOT_CINT_EXECUTABLE}
      ARGS -f ${DictNameCxx} -s ${TagLibName}
        ${INCLUDE_DIRS} ${IncludeFiles} ${LinkDefH}
      DEPENDS ${IncludeFiles} ${LinkDefH}
    )
  ENDIF()

ENDFUNCTION( ROOT6_GENERATE_DICTIONARY )
