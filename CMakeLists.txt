#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)
project(UCNsim)

#----------------------------------------------------------------------------
# Find Geant4 package, activating all available UI and Vis drivers by default
# You can set WITH_GEANT4_UIVIS to OFF via the command line or ccmake/cmake-gui
# to build a batch mode only executable
#
option(WITH_GEANT4_UIVIS "Build example with Geant4 UI and Vis drivers" ON)
if(WITH_GEANT4_UIVIS)
  find_package(Geant4 REQUIRED ui_all vis_all)
else()
  find_package(Geant4 REQUIRED)
endif()

#----------------------------------------------------------------------------
# Setup Geant4 include directories and compile definitions
#
include(${Geant4_USE_FILE})

# ASSIMP
find_package(assimp)
include_directories(/usr/local/assimp/include)
link_directories(/usr/local/assimp/lib)

# CADMesh
find_package(cadmesh)
include_directories(/usr/local/cadmesh/include)
link_directories(/usr/local/cadmesh/lib)

# muparser
include_directories(/usr/local/muparser/include)
link_directories(/usr/local/muparser/lib)

# alglib
file(GLOB alglib_sources ${PROJECT_SOURCE_DIR}/alglib-3.9.0/cpp/src/*.cpp)
file(GLOB alglib_headers ${PROJECT_SOURCE_DIR}/alglib-3.9.0/cpp/src/*.h)
include_directories(${PROJECT_SOURCE_DIR}/alglib-3.9.0/cpp/src)

# libtricubic
file(GLOB libtricubic_sources ${PROJECT_SOURCE_DIR}/libtricubic/*.cpp)
file(GLOB libtricubic_headers ${PROJECT_SOURCE_DIR}/libtricubic/*.h)
include_directories(${PROJECT_SOURCE_DIR}/libtricubic)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${Geant4_INCLUDE_DIR})
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the executable, and link it to the Geant4 libraries
#
#add_executable(UCNsim UCNsim.cc ${sources} ${headers})
add_executable(UCNsim UCNsim.cc ${sources} ${headers} ${alglib_sources} ${alglib_headers} ${libtricubic_sources} ${libtricubic_headers})
target_link_libraries(UCNsim ${Geant4_LIBRARIES})
target_link_libraries(UCNsim ${cadmesh_LIBRARIES})
target_link_libraries(UCNsim /usr/local/muparser/lib/libmuparser.so)

#----------------------------------------------------------------------------
# Install the executable to 'bin' directory under CMAKE_INSTALL_PREFIX
#
install(TARGETS UCNsim DESTINATION bin)
