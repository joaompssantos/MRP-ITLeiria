###
### CMakeLists for Minimum Rate Predictos
###

######################################################################################
# CMake  functions
######################################################################################
cmake_minimum_required(VERSION 3.10)

######################################################################################
# Project Definition
######################################################################################
set(PROJECT_NAME "4D-MRP")
project(${PROJECT_NAME})

set(MRP_VERSION_MAJOR "1")
set(MRP_VERSION_MINOR "1")
set(MRP_VERSION_PATCH "0")
set(MRP_VERSION "${MRP_VERSION_MAJOR}.${MRP_VERSION_MINOR}.${MRP_VERSION_PATCH}")
set(MRP_VERSION_STRING "${MRP_VERSION}")

set(MRP_VERSION_DATE "February 2019")

######################################################################################
# Variables
######################################################################################
OPTION(USE_EXTRA_WARNINGS "Enable extra warnings" OFF)
OPTION(USE_WERROR "Warnings as errors" OFF)

######################################################################################
# Configure targets
######################################################################################
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_BINARY_DIR})

# Define C standard 11
set(CMAKE_C_STANDARD 11)

# Define build type if not defined
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release")
    set(CMAKE_BUILD_TYPE "Release")
endif ()

# Set compile flags
set(CMAKE_C_FLAGS "-funroll-loops -Wall -Wextra" CACHE STRING "compile flags" FORCE)

# Set extra warnigns
if (USE_EXTRA_WARNINGS)
    string(CONCAT CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wfloat-equal -Wundef -Wshadow -Wpointer-arith -Wcast-align "
            "-Wstrict-prototypes -Wswitch-default -Wswitch-enum -Wconversion -Wunreachable-code")
endif ()

# Set warnings as errors
if (USE_WERROR)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror")
endif ()

# Add optimisation flags for DEBUG and RELEASE
set(CMAKE_C_FLAGS_DEBUG "-g -O0")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")

# Configure a header file to pass some of the CMake settings to the source code
configure_file("${PROJECT_SOURCE_DIR}/mrp_config.h.in" "${PROJECT_BINARY_DIR}/mrp_config.h")

set(COMMON_SRC
        mrp.h
        common.c
        entropy_coding.c
        rc.c)

# Add the executable
add_executable(encmrp.exe encmrp.c ${COMMON_SRC})
find_library(M_LIB m)
target_link_libraries(encmrp.exe ${M_LIB})

# Add the executable
add_executable(decmrp.exe decmrp.c ${COMMON_SRC})
find_library(M_LIB m)
target_link_libraries(decmrp.exe ${M_LIB})

######################################################################################
# Final information dump
######################################################################################
message("")
message(STATUS "${PROJECT_NAME} INFORMATION")
message(STATUS "Version: " "${MRP_VERSION_STRING}")
message(STATUS "Configuration:")
message(STATUS "    C++ Standard: " "${CMAKE_C_STANDARD}")
message(STATUS "    Build type: " "${CMAKE_BUILD_TYPE}")
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
    message(STATUS "    Build flags: " "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
else ()
    message(STATUS "    Build flags: " "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
endif ()
message("")