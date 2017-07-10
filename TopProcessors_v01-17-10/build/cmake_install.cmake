# Install script for directory: /home/aw/aw_Marlin/TopProcessors_v01-17-10

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/home/aw/aw_Marlin/TopProcessors_v01-17-10")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "RelWithDebInfo")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/." TYPE DIRECTORY FILES "/home/aw/aw_Marlin/TopProcessors_v01-17-10/./include" FILES_MATCHING REGEX "/[^/]*\\.h$" REGEX "/[^/]*\\~$" EXCLUDE REGEX "/[^/]*\\#[^/]*$" EXCLUDE REGEX "/\\.\\#[^/]*$" EXCLUDE REGEX "/[^/]*CVS$" EXCLUDE REGEX "/[^/]*\\.svn$" EXCLUDE)
endif()

if(NOT CMAKE_INSTALL_COMPONENT OR "${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified")
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libtopReconstruction_v17_10.so.0.1.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libtopReconstruction_v17_10.so.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libtopReconstruction_v17_10.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHECK
           FILE "${file}"
           RPATH "/home/aw/aw_Marlin/TopProcessors_v01-17-10/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/Marlin/v01-09/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/lcio/v02-07-02/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/gear/v01-06/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/CLHEP/2.1.4.1/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/ilcutil/v01-03/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/MarlinKinfit/v00-03/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/root/5.34.30/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/MarlinUtil/v01-12/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/CED/v01-09-01/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/gsl/2.1/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/RAIDA/v01-07/lib")
    endif()
  endforeach()
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE SHARED_LIBRARY PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE FILES
    "/home/aw/aw_Marlin/TopProcessors_v01-17-10/build/lib/libtopReconstruction_v17_10.so.0.1.0"
    "/home/aw/aw_Marlin/TopProcessors_v01-17-10/build/lib/libtopReconstruction_v17_10.so.0.1"
    "/home/aw/aw_Marlin/TopProcessors_v01-17-10/build/lib/libtopReconstruction_v17_10.so"
    )
  foreach(file
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libtopReconstruction_v17_10.so.0.1.0"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libtopReconstruction_v17_10.so.0.1"
      "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libtopReconstruction_v17_10.so"
      )
    if(EXISTS "${file}" AND
       NOT IS_SYMLINK "${file}")
      file(RPATH_CHANGE
           FILE "${file}"
           OLD_RPATH "/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/Marlin/v01-09/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/lcio/v02-07-02/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/gear/v01-06/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/CLHEP/2.1.4.1/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/ilcutil/v01-03/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/MarlinKinfit/v00-03/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/root/5.34.30/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/MarlinUtil/v01-12/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/CED/v01-09-01/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/gsl/2.1/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/RAIDA/v01-07/lib:::::::::::::::::::::::::::::::::::::::::::::::"
           NEW_RPATH "/home/aw/aw_Marlin/TopProcessors_v01-17-10/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/Marlin/v01-09/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/lcio/v02-07-02/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/gear/v01-06/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/CLHEP/2.1.4.1/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/ilcutil/v01-03/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/MarlinKinfit/v00-03/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/root/5.34.30/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/MarlinUtil/v01-12/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/CED/v01-09-01/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/gsl/2.1/lib:/cvmfs/ilc.desy.de/sw/x86_64_gcc44_sl6/v01-17-10/RAIDA/v01-07/lib")
      if(CMAKE_INSTALL_DO_STRIP)
        execute_process(COMMAND "/usr/bin/strip" "${file}")
      endif()
    endif()
  endforeach()
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/aw/aw_Marlin/TopProcessors_v01-17-10/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
