#cmake_minimum_required(VERSION 3.5)
cmake_minimum_required(VERSION 3.20.0)

project(qc2ga_sample)

# set the module path for "findLib.cmake"
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules/")

# require Eigen
find_package(Eigen3 REQUIRED)
if(${EIGEN3_FOUND})
    message(STATUS "lib EIGEN3 found")
    message(STATUS "  version " ${EIGEN3_VERSION_STRING})
    message(STATUS "  include " ${EIGEN3_INCLUDE_DIR})
endif()


# require lib qc2ga
find_package(QC2GA)
if(${QC2GA_FOUND})
    message(STATUS "lib QC2GA found")
    message(STATUS "  include " ${QC2GA_INCLUDE_DIRS})
    message(STATUS "  lib     " ${QC2GA_LIBRARIES})
else()
    message(STATUS "lib QC2GA should be set manually, please uncomment the commented lines in the CMakeLists.txt file")
    #set(QC2GA_INCLUDE_DIRS "[path to]/garamon_qc2ga/src/")
    #set(QC2GA_LIBRARIES "[path to]/garamon_qc2ga/build/libqc2ga.dylib")

    # for manual setting, check if valid
    if(EXISTS "${QC2GA_INCLUDE_DIRS}")
        message(STATUS "include for lib qc2ga found: ${QC2GA_INCLUDE_DIRS}" )
    else()
        message(FATAL_ERROR "include for lib qc2ga invalid or not found: ${QC2GA_INCLUDE_DIRS}" )
    endif()

    # for manual setting, check if valid
    if(EXISTS "${QC2GA_LIBRARIES}")
        message(STATUS "library qc2ga found: ${QC2GA_LIBRARIES}" )
    else()
        message(FATAL_ERROR "library qc2ga not found: ${QC2GA_LIBRARIES}" )
    endif()
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -std=c++20")

add_executable(solve_inter src/solve_conic_intersection.cpp src/conic.hpp )

# includes
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${QC2GA_INCLUDE_DIRS})


target_link_libraries(solve_inter ${QC2GA_LIBRARIES})
