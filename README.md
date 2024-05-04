---------------------------------------------------------------------------
## IMPORTANT

Before running cmake
1. [CERN ROOT](https://root.cern.ch) must be installed.\
   The `root-config` script must be in `$ROOTSYS/bin` directory,
   otherwise, you must specify the path to `root-config` when call
   the cmake.
   For example, in IHEP cluster (CentOS7):
   ```
   % cmake -DROOT_CONFIG_SEARCHPATH="/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/lcg/LCG_84/ROOT/6.20.02/x86_64-centos7-gcc49-opt/bin"
   ```

2. [CLHEP](https://proj-clhep.web.cern.ch/proj-clhep) library must be
   installed.\
   The path to CLHEP can be specified in the environment variable
   `CLHEP_DIR`. Also you can specify the path to CLHEP when calling
   cmake.
   For example, in IHEP cluster (CentOS7):
   ```
   % cmake -DCLHEP_SEARCHPATH="/cvmfs/bes3.ihep.ac.cn/bes3sw/ExternalLib/lcg/LCG_84/clhep/2.4.4.0/x86_64-centos7-gcc49-opt"
   ```

---------------------------------------------------------------------------
## INSTALLATION

1. Create build directory: ```% mkdir build; cd build```

2. Run `cmake ../ [cmake-options]` to prepare the infrastructure
   for building the BEAN program.
   * For example:\
     `% cmake ../ -DBOSS_VERSION=6.6.4`\
     or\
     `% cmake ../ -DBOSS_VERSION=7.0.4`
   * If you want to compile program using a specific compiler,
     use the following command:\
     `% cmake -DCMAKE_C_COMPILER=gcc49 -DCMAKE_CXX_COMPILER=g++49 ../ [cmake-options]`
   * If you wish, you can run `ccmake` and change the options
     in the dialog interface.

3. Compile the BEAN by runing `% make` and then `% make install`
   * To see detailed output of compilation commands, use
     `% make VERBOSE=1`
   * `make install` will install program, libraries and script
     in the `workdir/` directory

4. The BEAN program uses sqlite database located in the
   `Analysis/DatabaseSvc/dat/` directory. Currently, updating this
   database is only possible "manually". You need to go to the
   `http://docbes3.ihep.ac.cn/db_ana/` website, download and unzip
   three files: `db.timestamp  offlinedb.db  run.db`.\
   _Note: The `db.timestamp` file contains Unix time and md5 checksums
   for database files with which you can verify the integrity of
   downloaded files._

5. To remove BEAN files from the `workdir/` directory, use the `% make
   uninstall` command. This may be useful when updating the program.
   Alternatively, you can delete the 'bean_{BOSS_VERSION}.exe' file
   and the 'BeanLib_{BOSS_VERSION}' directory manually.

6. After running make install, you can remove the build directory
   at any time.

---------------------------------------------------------------------------
## NOTES AND NEWS

* The user source programs MUST be in directory `BeanUser/`.\
  Edit `BeanUser/CMakeLists.txt` to add name of your file
  to the list of files to be compiled.

* The AbsCor algorithm in BEAN does not use the database, unlike
  the BOSS version. However, you can read calibration files using
  the functions:
  ```
    void ReadDatac3p(std::string DataPathc3p, bool info=true);
    void ReadParMcCor(std::string paraPath, bool info=true);
    void ReadCorFunpara(std::string CorFunparaPath, bool info=true);
  ```
  The path to the calibration files in CernVM `/cvmfs/` file system
  can be obtained using the python3 program `GetAbsCorFiles.py` located
  in `bean/scripts` directory.\
  An example of using these functions is contained in the file:
  `BeanUser/TestAbsCor.cxx`

* A small bug in the VertexFit package where the total energy was not
  updated after fitting has been fixed for all versions of the BEAN
  program.

* The `DecayTable` class has been added to `EventTag` algorithm.
  This class is designed to fill out a Monte Carlo particle decay
  table in text form.

* The `RscanDQ` package has been added to the BEAN algorithms.
  When you using the 2015 Rscan data, you must use this package
  to exclude bad runs. For more information, see the
  Tau and QCD Group Data Samples page:
  (https://docbes3.ihep.ac.cn/~tauqcdgroup/index.php/Data_Samples)

* The `TrackCorrection` class for the helix parameters corrections
  for MC tracks, based on "TrackCorrection" package has been added
  to the BEAN algorithms.
  For more information, see the Charmonium Group's Data Quality page:
  (https://docbes3.ihep.ac.cn/~charmoniumgroup/index.php/DataQuality_Page)

---------------------------------------------------------------------------
## DOCUMENTATION

* Somewhat outdated, but still useful documentation on the BEAN
  on the BES3 Offline Software Group page:
  (https://docbes3.ihep.ac.cn/~offlinesoftware/index.php/BEAN)

* Examples in the BeanUser directory:

| File             | Description                                       |
| :---             | :---                                              |
| ---------------- | ------------------------------------------------- |
| UserTest.cxx     | two simple examples of writing user functions in  |
| User1.cxx        | the BEAN; demonstration of saving histograms      |
| ---------------- | ------------------------------------------------- |
| Rhopi.cxx        | Rhopi program from `Analysis/Physics/RhopiAlg`    |
| ---------------- | ------------------------------------------------- |
| TestAbsCor.cxx   | example of using functions of AbsCor algorithms   |
| TestPID.cxx      | example of using the ParticleID library functions |
| MagField.cxx     | example of using the MagneticField functions      |
| TestEventTag.cxx | example of using the EventTag facility            |
| TestDb.cxx       | example of using the SqliteInterface class        |
| ---------------- | ------------------------------------------------- |
| JobInfo.cxx      | program for printing information from JobInfoTree |
| EntryList.cxx    | program with an example of using TEntryList       |
| ---------------- | ------------------------------------------------- |

---------------------------------------------------------------------------
## TODO

* Add an example of using the `DecayTable` class in `TestEventTag.cxx`

* The installation and operation of BEAN under Windows has been tested
  for a very long time. I have not tested it at all with the version of
  ROOT-6 under Windows.

* `% make proofbean` command creates PAR files for ROOT-PROOF.
  But the work of BEAN on the PROOF cluster has not been tested for a very
  long time.

* `% make setup_file` command creates the `setup_{BOSS_VERSION}.sh` file.
  This bash script sets environment variables with paths to libraries used
  in the BEAN. In the case of a standard installation of all components,
  this script **IS NOT REQUIRED**.
  Perhaps it could be useful for debugging.

* `% make updatedb` command was intended to be used to install or update
  a database. After entering password access, this method no longer (yet?)
  works.

---------------------------------------------------------------------------
