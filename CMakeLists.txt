cmake_minimum_required(VERSION 3.5.1 FATAL_ERROR)
#set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
#set(CMAKE_MACOSX_RPATH TRUE)

macro(DATE RESULT)
  execute_process(COMMAND "date" "+%Y-%m-%d" OUTPUT_VARIABLE ${RESULT}
                                             OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()

macro(TIME RESULT)
  execute_process(COMMAND "date" "+%H:%M:%S" OUTPUT_VARIABLE ${RESULT}
                                             OUTPUT_STRIP_TRAILING_WHITESPACE)
endmacro()

# MSHMINT EXECUTABLE
project(mshint LANGUAGES C)
set(CMAKE_BUILD_TYPE "Release")

# Write compilation date in a compil.date file of the sources/mshint directory
DATE(day)
TIME(hour)
file(WRITE sources/commonFiles/compil.date "#define COMPIL \"${day} ${hour}\"")

# Add the source files to the mshint executable
add_executable(mshint sources/mshint/boule.c
                      sources/mshint/bucket.c
                      sources/commonFiles/chrono.c
                      sources/mshint/hash.c
                      sources/mshint/inout.c
                      sources/mshint/intelt.c
                      sources/commonFiles/libmesh5.c
                      sources/mshint/locelt.c
                      sources/commonFiles/memory.c
                      sources/mshint/mshin1.c
                      sources/mshint/mshint.c
                      sources/mshint/scalem.c)

# Find the math library
find_library(M_LIB m)

# Set all targets and properties of the mshint executable
target_link_libraries(mshint PUBLIC ${M_LIB})
target_include_directories(mshint PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources/mshint>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources/commonFiles>)
set_property(TARGET mshint PROPERTY C_STANDARD 99)

# Set the installation properties of the mshint executable
install(TARGETS mshint RUNTIME DESTINATION "${PROJECT_SOURCE_DIR}/bin")
install(FILES sources/commonFiles/chrono.h
              sources/commonFiles/libmesh5.h
              sources/commonFiles/memory.h
              sources/mshint/mshint.h DESTINATION "${PROJECT_SOURCE_DIR}/include")


# MSHMET EXECUTABLE
project(mshmet LANGUAGES C)
set(CMAKE_BUILD_TYPE "Release")

# Add the source files to the mshmet executable
add_executable(mshmet sources/mshmet/boulep.c
                      sources/mshmet/defmet.c
                      sources/commonFiles/eigenv.c
                      sources/mshmet/gradls.c
                      sources/commonFiles/chrono.c
                      sources/mshmet/hash.c
                      sources/mshmet/hessls.c
                      sources/mshmet/inout.c
                      sources/mshmet/laplac.c
                      sources/commonFiles/libmesh5.c
                      sources/commonFiles/memory.c
                      sources/mshmet/mshme1.c
                      sources/mshmet/mshmet.c
                      sources/mshmet/redsim.c
                      sources/mshmet/scalem.c)

# Find the math library
find_library(M_LIB m)

# Set all targets and properties of the mshmet executable
target_link_libraries(mshmet PUBLIC ${M_LIB})
target_include_directories(mshmet PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources/mshmet>
                                         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sources/commonFiles>)
set_property(TARGET mshmet PROPERTY C_STANDARD 99)

# Set the installation properties of the mshmet executable (only the non-already copied files)
install(TARGETS mshmet RUNTIME DESTINATION "${PROJECT_SOURCE_DIR}/bin")
install(FILES #sources/commonFiles/chrono.h
              sources/commonFiles/eigenv.h
              #sources/commonFiles/libmesh5.h
              #sources/commonFiles/memory.h
              sources/mshmet/mshmet.h DESTINATION "${PROJECT_SOURCE_DIR}/include")

