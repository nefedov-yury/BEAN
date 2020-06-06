------------------------------------------------------------------------------
 ##IMPORTANT

 Before running cmake
 A) CERN ROOT must be installed
        https://root.cern.ch
    'root-config' must be in '$ROOTSYS/bin/' directory,
    otherwise, you must specify the path to 'root-config'
    when call the cmake:

        cmake -DROOT_CONFIG_SEARCHPATH="/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/external/ROOT/5.34.09/x86_64-slc6-gcc46-opt/root/bin"

 B) CLHEP library must be installed
        https://proj-clhep.web.cern.ch/proj-clhep/
    The path to CLHEP may be specified in the environment
    variable CLHEP_DIR or when call the cmake:

        cmake -DCLHEP_SEARCHPATH="/afs/ihep.ac.cn/bes3/offline/ExternalLib/SLC6/ExternalLib/external/clhep/2.0.4.5/x86_64-slc6-gcc46-opt"

------------------------------------------------------------------------------
 ##INSTALLATION

 1) create build directory (mkdir build; cd build)

 2) cmake ../ ( + optionally ccmake ../ )
    Examples:
       cmake ../ -DBOSS_VERSION=6.6.4
    or
       cmake ../ -DBOSS_VERSION=7.0.4

    if you want to compile with a specific version of the compiler:
       cmake -DCMAKE_C_COMPILER=gcc46 -DCMAKE_CXX_COMPILER=g++46 ../

 3) make
        to see compilation commands use 'make VERBOSE=1'
    make install
        this will install program, libraries and script in the
        workdir/ directory
    make uninstall
        this will remove program files from workdir/;
        this may be useful when updating the program

 4) to install or update the databases use the following command
    make updatedb
        It requires python and wget. In case of their absence,
        you have to install the databases in some other way.
        For example, you may copy them from another machine.
        The path to database is Analysis/DatabaseSvc/dat/

 5) to create PAR-files for ROOT-PROOF use
    make proofbean

------------------------------------------------------------------------------
  ###NOTE

  The executable file "bean_{BOSS_VERSION}.exe"
  is installed in the "workdir" directory.
  After installing the program, you can delete the build
  directory at any time, all libraries and executable are
  in the workdir/ directory.

  The user source programs MUST be in directory "BeanUser":
  edit BeanUser/CMakeLists.txt to add name of your file
  to the list of files to be compiled.

  Script "setup_{BOSS_VERSION}.sh" contains settings
  for the environment variable LD_LIBRARY_PATH;
  In the case of a normal installation of the program,
  this script IS NOT REQUIRED.

------------------------------------------------------------------------------
  ###NEWS

  1. RscanDQ: When you using the 2015 Rscan data, you must use this package
     to exclude bad runs:
     (https://docbes3.ihep.ac.cn/~tauqcdgroup/index.php/Data_Samples)

  2. New functionality in VertexDbSvc: new function
     ReadOneTime(int runFrom, int runTo) - you can get verteces from
     database at once for runs in the interval [runFrom,runTo].
     See BeanUser/TestDb.cxx how to use this function

------------------------------------------------------------------------------
