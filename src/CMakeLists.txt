cmake_minimum_required(VERSION 3.8)
project(src)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES
        flowpairgenerator.hpp
        helpers.hpp
        ksp.hpp
        main.cpp
        testcases.hpp
        two_flows_update.hpp)

add_executable(src ${SOURCE_FILES})