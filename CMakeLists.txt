cmake_minimum_required(VERSION 3.6)
project(FLEX_MRP)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(ENC_SOURCE_FILES
        common.c
        encmrp.c
        mrp.h
        rc.c)

add_executable(encmrp.exe ${ENC_SOURCE_FILES})
find_library(M_LIB m)
target_link_libraries(encmrp.exe ${M_LIB})

set(DEC_SOURCE_FILES
        common.c
        decmrp.c
        mrp.h
        rc.c)

add_executable(decmrp.exe ${DEC_SOURCE_FILES})
find_library(M_LIB m)
target_link_libraries(decmrp.exe ${M_LIB})