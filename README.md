---------------------------------------------------------------------------
## IMPORTANT
Before running cmake
1. [CERN ROOT](https://root.cern.ch) must be installed.\
   The `root-config` script must be in `$ROOTSYS/bin` directory, otherwise,
   you must specify the path to `root-config` when call the cmake.
   For example, in IHEP cluster (CentOS7):
```
% cmake -DROOT_CONFIG_SEARCHPATH="/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/lcg/LCG_84/ROOT/6.20.02/x86_64-centos7-gcc49-opt/bin"
```

2. [CLHEP](https://proj-clhep.web.cern.ch/proj-clhep) library must be
   installed.\
   The path to CLHEP can be specified in the environment variable
   `CLHEP_DIR`. Also you can specify the path to CLHEP when calling cmake.
   For example, in IHEP cluster (CentOS7):
```
% cmake -DCLHEP_SEARCHPATH="/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/lcg/LCG_84/clhep/2.4.4.0/x86_64-centos7-gcc49-opt"
```

---------------------------------------------------------------------------
## INSTALLATION
1. Create build directory: ```% mkdir build; cd build```

2. Run `cmake ../ [cmake-options]` to prepare the infrastructure
   for building the BEAN program
   * For example:\
     `% cmake ../ -DBOSS_VERSION=6.6.4`\
     or\
     `% cmake ../ -DBOSS_VERSION=7.0.4`
   * If you want to compile program with a specific compiler version,
     use the following command:\
     `% cmake -DCMAKE_C_COMPILER=gcc49 -DCMAKE_CXX_COMPILER=g++49 ../ [cmake-options]`
   * If you wish, you can run ``ccmake`` and change the options in the
     dialog interface.

3. Compile BEAN by runing `% make` and then `% make install`
   * To see detailed output of compilation commands, use
     `% make VERBOSE=1`
   * `make install` will install program, libraries and script in the
     `workdir/` directory
   * To remove BEAN files from the `workdir/` directory,
     use the `% make uninstall` command. This may be useful when
     updating the program.

4. Use the following command to install or update the database:
   `% make updatedb`\
    _It requires python and wget. In case of their absence,
    you have to install the databases in some other way.
    For example, you may copy them from another machine.
    The path to database is_ `Analysis/DatabaseSvc/dat/`

5. PAR files for ROOT-PROOF can be created with the
   `% make proofbean` command.\
   _**NOTE:**_ This option has not been tested for a very long time

---------------------------------------------------------------------------
## NOTE

* After installing BEAN, you can delete the build directory at any time.

* The user source programs MUST be in directory `BeanUser/`.
  Edit `BeanUser/CMakeLists.txt` to add name of your file
  to the list of files to be compiled.

* The `setup_{BOSS_VERSION}.sh` script sets environment variables
   with paths to libraries used in BEAN. In the case of a standard
   installation of all components, this script **IS NOT REQUIRED**.

---------------------------------------------------------------------------
## NEWS

1. 'RscanDQ': When you using the 2015 Rscan data, you must use this
   package to exclude bad runs. For details see:
```
   https://docbes3.ihep.ac.cn/~tauqcdgroup/index.php/Data_Samples
```

2. The `TrackCorrection` class based on "TrackCorrection" package
   modified for use in a BEAN (see `~guoyp/public/TrackCorrection`).
   For details, please refer the Data Qality page:
```
   https://docbes3.ihep.ac.cn/~charmoniumgroup/index.php/DataQuality_Page
```

------------------------------------------------------------------------------
