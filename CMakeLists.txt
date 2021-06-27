cmake_minimum_required (VERSION 3.6)

SET(EXE_NAME "WFCTiles")

project (${EXE_NAME})

file(GLOB_RECURSE SRC
    src/*.h
    src/*.cpp
)



ADD_EXECUTABLE(${EXE_NAME} ${SRC})

source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR} FILES ${SRC})