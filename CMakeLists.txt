cmake_minimum_required(VERSION 2.6)
project(JpegDecoder)

set(SRC_DIR "${PROJECT_SOURCE_DIR}/src")
include_directories(${SRC_DIR})
add_subdirectory(src)

file(GLOB JPEG_DECODER_SRC
     ${SRC_DIR}/*.h
     ${SRC_DIR}/*.cpp
     )
add_executable(JpegDecoder main.cpp ${JPEG_DECODER_SRC})