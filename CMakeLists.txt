cmake_minimum_required(VERSION 3.6)
project(MRP)

#set(CMAKE_BUILD_TYPE "Release")

IF(UNIX)
    IF(CMAKE_COMPILER_IS_GNUCC)
        IF(CMAKE_BUILD_TYPE MATCHES Debug)
            message("Debug build.")
            set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -g -Wall -funroll-loops")
            set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0 -g -Wall -funroll-loops")
        ELSEIF(CMAKE_BUILD_TYPE MATCHES Release)
            message("Release build.")
            set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -Wall -funroll-loops")
            set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} -O3 -Wall -funroll-loops")
        ENDIF()
    ENDIF(CMAKE_COMPILER_IS_GNUCC)
ENDIF(UNIX)

# Set binary path
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../)

# The version number.
#set (MRP_VERSION_MAJOR 0)
#set (MRP_VERSION_MINOR 9)

# Cnfigure a header file to pass some of the CMake settings
# to the source code
#configure_file (
#    "${PROJECT_SOURCE_DIR}/mrp.h.in"
#    "${PROJECT_SOURCE_DIR}/mrp.h"
#)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(ENC_SOURCE_FILES
        common.c
        encmrp.c
        mrp.h
        rc.c)

# Add the executable
add_executable(encmrp.exe ${ENC_SOURCE_FILES})
find_library(M_LIB m)
target_link_libraries(encmrp.exe ${M_LIB})

set(DEC_SOURCE_FILES
        common.c
        decmrp.c
        mrp.h
        rc.c)

# Add the executable
add_executable(decmrp.exe ${DEC_SOURCE_FILES})
find_library(M_LIB m)
target_link_libraries(decmrp.exe ${M_LIB})