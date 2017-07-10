########################################################
# cmake file for building Marlin example Package
# @author Jan Engels, Desy IT
CMAKE_MINIMUM_REQUIRED(VERSION 2.6 FATAL_ERROR)
########################################################


# project name
PROJECT( topReconstruction_v17_10 )


# project version
SET( ${PROJECT_NAME}_VERSION_MAJOR 0 )
SET( ${PROJECT_NAME}_VERSION_MINOR 1 )
SET( ${PROJECT_NAME}_VERSION_PATCH 0 )



### DEPENDENCIES ############################################################

FIND_PACKAGE( ILCUTIL REQUIRED COMPONENTS ILCSOFT_CMAKE_MODULES )

# load default settings from ILCSOFT_CMAKE_MODULES
INCLUDE( ilcsoft_default_settings )


FIND_PACKAGE( Marlin 1.0 REQUIRED ) # minimum required Marlin version
#FIND_PACKAGE( MarlinUtil 1.4 REQUIRED )
FIND_PACKAGE( MarlinReco )
FIND_PACKAGE( MarlinKinfit )
FIND_PACKAGE( ROOT )
#aw 15/1/16 added MarlinUtil to accomodate TauFinder Processor
FIND_PACKAGE( MarlinUtil )

FOREACH( pkg Marlin MarlinReco MarlinKinfit ROOT MarlinUtil)
  IF( ${pkg}_FOUND )
      INCLUDE_DIRECTORIES( ${${pkg}_INCLUDE_DIRS} )
      LINK_LIBRARIES( ${${pkg}_LIBRARIES} )
      ADD_DEFINITIONS( ${${pkg}_DEFINITIONS} )
  ENDIF()
ENDFOREACH()


FIND_PACKAGE( CLHEP REQUIRED )
INCLUDE_DIRECTORIES( ${CLHEP_INCLUDE_DIRS} )
LINK_LIBRARIES( ${CLHEP_LIBRARIES} )

FIND_PACKAGE( GSL 1.12 REQUIRED )
INCLUDE_DIRECTORIES( ${GSL_INCLUDE_DIRS} )
LINK_LIBRARIES( ${GSL_LIBRARIES} )

# optional package
FIND_PACKAGE( AIDA )
IF( AIDA_FOUND )
    INCLUDE_DIRECTORIES( ${AIDA_INCLUDE_DIRS} )
    LINK_LIBRARIES( ${AIDA_LIBRARIES} )
    ADD_DEFINITIONS( "-DMARLIN_USE_AIDA" )
    MESSAGE( STATUS "AIDA -- found" )
ELSE()
    MESSAGE( STATUS "AIDA -- not found" )
ENDIF()




### DOCUMENTATION ###########################################################

OPTION( INSTALL_DOC "Set to OFF to skip build/install Documentation" OFF )

IF( INSTALL_DOC AND EXISTS "${PROJECT_SOURCE_DIR}/doc/CMakeLists.txt" )
    ADD_SUBDIRECTORY( ./doc )
ENDIF()



### LIBRARY #################################################################

# definitions to pass to the compiler
#ADD_DEFINITIONS( "-Wall -ansi -std=c++11" )
ADD_DEFINITIONS( "-Wall -ansi " )
ADD_DEFINITIONS( "-Wno-long-long" )

# include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL_DIRECTORY( ./include DESTINATION . FILES_MATCHING PATTERN "*.h" )

# add library
AUX_SOURCE_DIRECTORY( ./src library_sources )
ADD_SHARED_LIBRARY( ${PROJECT_NAME} ${library_sources} )
INSTALL_SHARED_LIBRARY( ${PROJECT_NAME} DESTINATION lib )


# display some variables and write them to cache
DISPLAY_STD_VARIABLES()



