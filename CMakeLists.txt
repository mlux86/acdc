cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(unnamed)
enable_testing()

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -std=c++11")

# Boost
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.55.0 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# SeqAn
set(SEQAN_FIND_DEPENDENCIES ZLIB BZip2)
find_package(SeqAn REQUIRED)
include_directories(${SEQAN_INCLUDE_DIRS})
add_definitions(${SEQAN_DEFINITIONS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

add_subdirectory(libdip)
add_subdirectory(libtsne)
add_subdirectory(src)
add_subdirectory(test)
