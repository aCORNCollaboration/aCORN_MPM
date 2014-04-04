#----------------------------------------------------------------------------
# Usage:
# cmake -DGeant4_DIR=~/geant4.9.6/geant4.9.6-install/lib/Geant4-9.6.2/ ../aCORN_G4/
# make -j6

#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(ucnG4_prod)

# build libaCORN_MPM.a dependency
add_custom_target(
	extern_lib
	COMMAND make
	WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/..
)


#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)
include_directories(${PROJECT_SOURCE_DIR}/include
					${PROJECT_SOURCE_DIR}/../Physics
					${PROJECT_SOURCE_DIR}/../IOUtils
					${PROJECT_SOURCE_DIR}/../RootUtils)
link_directories(${PROJECT_SOURCE_DIR}/..)

#-------------------------------
# ROOT dependencies
#-------------------------------

######## locate ROOT install, headers, includes
include(FindROOT.cmake)
include_directories(${ROOT_INCLUDE_DIR})
link_directories(${ROOT_LIBRARY_DIR})

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()
include(${Geant4_USE_FILE})


##############
# Output paths
##############

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")


###############
# TObject subclasses dictionary & aCORN_G4 shared library
###############

set(aCORN_G4_LINKDEF ${PROJECT_SOURCE_DIR}/aCORN_G4_LinkDef.h)
set(aCORN_G4_DICTIONARY ${CMAKE_CURRENT_BINARY_DIR}/aCORN_G4_Dict.cxx)
MESSAGE("Generating dictionary ${aCORN_G4_DICTIONARY}")
set(DICTFILES "${PROJECT_SOURCE_DIR}/include/MCEvent.hh;${PROJECT_SOURCE_DIR}/include/TrackInfo.hh;${PROJECT_SOURCE_DIR}/include/PrimaryInfo.hh")
ROOT_GENERATE_DICTIONARY("${DICTFILES}" "${aCORN_G4_LINKDEF}" "${aCORN_G4_DICTIONARY}" "")
LIST(APPEND sources ${aCORN_G4_DICTIONARY})

add_library(aCORN_G4 SHARED ${sources} ${headers})
add_dependencies(aCORN_G4 extern_lib)
target_link_libraries(aCORN_G4 ${ROOT_LIBRARIES} ${Geant4_LIBRARIES} aCORN_MPM)
SET(ucnG4_LIBRARY_PROPERTIES
		${aCORN_G4_LIBRARY_PROPERTIES}
		VERSION "0"
		SOVERSION "0"
		SUFFIX ".so"
)
set_target_properties(ucnG4 PROPERTIES ${ucnG4_LIBRARY_PROPERTIES})
install(TARGETS ucnG4 DESTINATION ${LIBRARY_OUTPUT_PATH})
link_directories(${LIBRARY_OUTPUT_PATH})


#################
# Executables
#################

add_executable(aCORN_G4_prod aCORN_G4_prod.cc)
target_link_libraries(aCORN_G4_prod aCORN_G4 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} libaCORN_MPM.a)
add_dependencies(ucnG4_prod extern_lib aCORN_G4)

add_executable(aCORN_MC_Analyzer aCORN_MC_Analyzer.cc)
target_link_libraries(aCORN_MC_Analyzer aCORN_G4 ${Geant4_LIBRARIES} ${ROOT_LIBRARIES} libaCORN_MPM.a)
add_dependencies(aCORN_MC_Analyzer extern_lib aCORN_G4)

install(TARGETS aCORN_G4_prod aCORN_MC_Analyzer DESTINATION ${EXECUTABLE_OUTPUT_PATH})


#################
# Build settings
#################

#set(CMAKE_BUILD_TYPE debug)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
set(CMAKE_BUILD_TYPE release)

