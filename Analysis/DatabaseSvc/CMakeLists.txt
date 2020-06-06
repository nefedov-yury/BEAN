#------------------------------------------------------------------------------
# DatabaseSvc: interface library to sqlite database
#              <=> boss/Database/DatabaseSvc/
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# debug
#------------------------------------------------------------------------------
# PROJECT( db )
# CMAKE_MINIMUM_REQUIRED( VERSION 2.6.0 FATAL_ERROR )
# OPTION( BUILTIN_SQLITE "Use included sqlite, or system sqlite library" ON )
# ADD_DEFINITIONS( -DBEAN )
#------------------------------------------------------------------------------

SET( libname DatabaseSvc )

INCLUDE_DIRECTORIES( ${CMAKE_CURRENT_SOURCE_DIR} )

#------------------------------------------------------------------------------
# the list of source files
FILE( GLOB CXX_FILES src/*.cxx )

#------------------------------------------------------------------------------
# sqlite library

UNSET( SQLIB )
IF( BUILTIN_SQLITE )

  IF( NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/sqlite3.c OR
      NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/sqlite3.h )
      MESSAGE( FATAL_ERROR
         "There are no local files \"sqlite3.c\" and \"sqlite3.h\".\n"
         "Please download them from site: \nhttp://www.sqlite.org \n"
         "and put them in ${CMAKE_CURRENT_SOURCE_DIR} directory "
         "(see script Analysis/DatabaseSvc/share/get_sqlite3.sh).\n"
         "Also you can try to use the system sqlite library "
         "(set BUILTIN_SQLITE to OFF and rerun cmake)."
      )
  ENDIF()

  # add sqlite3.c to the list of source files
  SET( CXX_FILES ${CXX_FILES} ${CMAKE_CURRENT_SOURCE_DIR}/sqlite3.c )
#   MESSAGE( " +++DatabaseSvc: CXX_FILES= ${CXX_FILES}" )

  IF ( CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR
       CMAKE_CXX_COMPILER_ID MATCHES "GNU"       )
     SET_SOURCE_FILES_PROPERTIES( ${CMAKE_CURRENT_SOURCE_DIR}/sqlite3.c
       PROPERTIES COMPILE_FLAGS -Wno-deprecated-declarations )
  ENDIF()

  # There are various methods for file locking:
  #... 2. No locking (only for  Mac OS X)
  # ADD_DEFINITIONS( -DSQLITE_ENABLE_LOCKING_STYLE=2 )

ELSE()

  # try to find SQLITE3 library
  FIND_LIBRARY( SQLIB sqlite3 )
  IF( NOT SQLIB )
      MESSAGE( FATAL_ERROR
         "The system sqlite library was not found. "
         "Try to set BUILTIN_SQLITE to ON and rerun cmake."
      )
  ENDIF()

ENDIF( BUILTIN_SQLITE )

#------------------------------------------------------------------------------
# updatedb
IF( ${CMAKE_VERSION} VERSION_LESS "3.12.0" )
  find_package( PythonInterp )
  IF( PYTHONINTERP_FOUND )
    ADD_CUSTOM_TARGET( updatedb
      COMMAND "${CMAKE_COMMAND}" -E make_directory "../dat"
      COMMAND "${PYTHON_EXECUTABLE}" "./dbupdate2.py" "../dat"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/share/"
    )
  ENDIF()
ELSE()  
  find_package (Python3 COMPONENTS Interpreter)
  IF( Python3_Interpreter_FOUND )
    ADD_CUSTOM_TARGET( updatedb
      COMMAND "${CMAKE_COMMAND}" -E make_directory "../dat"
      COMMAND "${Python3_EXECUTABLE}" "./dbupdate.py" "../dat"
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/share/"
    )
  ELSE()
    IF( Python2_Interpreter_FOUND )
      ADD_CUSTOM_TARGET( updatedb
        COMMAND "${CMAKE_COMMAND}" -E make_directory "../dat"
        COMMAND "${Python2_EXECUTABLE}" "./dbupdate2.py" "../dat"
        WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/share/"
      )
    ENDIF()
  ENDIF()
ENDIF()
#------------------------------------------------------------------------------
#
IF( CMAKE_SYSTEM_NAME MATCHES Windows ) # WINDOWS
  ADD_LIBRARY( ${libname} STATIC ${CXX_FILES} )
ELSE()
  ADD_LIBRARY( ${libname} SHARED ${CXX_FILES} )
  MAC_LIB_RPATH( ${libname} )
ENDIF()

# We can not assume that pthread library will be included
# automatically with "root-config"
FIND_PACKAGE( Threads )
IF( SQLIB )
  TARGET_LINK_LIBRARIES( ${libname} ${SQLIB}
                                    ${CMAKE_THREAD_LIBS_INIT} )
ELSE()
  TARGET_LINK_LIBRARIES( ${libname} ${CMAKE_THREAD_LIBS_INIT} )
ENDIF()

#------------------------------------------------------------------------------
# install
#------------------------------------------------------------------------------
INSTALL( TARGETS ${libname} DESTINATION ${install_lib_dir} )

