cmake_minimum_required(VERSION 2.8)
set(CMAKE_BUILD_TYPE  "Release")
INCLUDE_DIRECTORIES(  "$ENV{HOME}/include")
LINK_DIRECTORIES(     "$ENV{HOME}/lib")
SET(CMAKE_CXX_FLAGS   "-Wunused-result -Wno-dev -g -O3")
set(CMAKE_MACOSX_RPATH 1)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#MshInt
project(_mshint)
#library
file(                  GLOB_RECURSE source_files mshint/*)
file(                  GLOB_RECURSE header_files mshint/*.h)
ADD_LIBRARY(           MshInt SHARED ${source_files})
find_library(          Commons NAMES Commons HINTS "$ENV{HOME}/lib")
target_link_libraries( MshInt ${Commons})
INSTALL(               TARGETS MshInt LIBRARY  DESTINATION "$ENV{HOME}/lib")
INSTALL(               FILES ${header_files} DESTINATION "$ENV{HOME}/include")
#Executable
add_executable(        mshint mshint/mshint.c)
target_link_libraries( mshint Commons MshInt)
INSTALL(               TARGETS mshint RUNTIME DESTINATION "$ENV{HOME}/bin")
set_property(TARGET    mshint PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)

#MshMet
project(_mshmet)
#library
file(                  GLOB_RECURSE source_files mshmet/*)
file(                  GLOB_RECURSE header_files mshmet/*.h)
ADD_LIBRARY(           MshMet SHARED ${source_files})
find_library(          Commons NAMES Commons HINTS "$ENV{HOME}/lib")
target_link_libraries( MshMet ${Commons})
INSTALL(               TARGETS MshMet LIBRARY  DESTINATION "$ENV{HOME}/lib")
INSTALL(               FILES ${header_files} DESTINATION "$ENV{HOME}/include")
#Executable
add_executable(        mshmet mshmet/mshmet.c)
target_link_libraries( mshmet Commons MshMet)
INSTALL(               TARGETS mshmet RUNTIME DESTINATION "$ENV{HOME}/bin")
set_property(TARGET    mshmet PROPERTY INSTALL_RPATH_USE_LINK_PATH TRUE)
