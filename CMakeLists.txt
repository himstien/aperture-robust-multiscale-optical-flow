# Copyright: 2018
# Author: Himanshu Akolkar


cmake_minimum_required(VERSION 2.6)
set (CMAKE_CXX_STANDARD 11)

set(PROJECTNAME FARMS_Flow)
project(${PROJECTNAME})


file(GLOB source "src/*.cpp")
file(GLOB header "include/*.h")

source_group("Source Files" FILES ${source})
source_group("Header Files" FILES ${header})

find_package (Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIRS})

message(STATUS "print")

add_executable(${PROJECTNAME} ${source} ${header})

target_link_libraries (${PROJECTNAME} ${EIGEN3_LIBRARIES})

install(TARGETS ${PROJECTNAME} DESTINATION bin)
