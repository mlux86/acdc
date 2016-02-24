cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(acdc)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -std=c++11")

# Boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.58.0 COMPONENTS program_options system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# Eigen
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIRS})

# t-SNE
include_directories(SYSTEM "libtsne")

# Dip statistic
include_directories(SYSTEM "libdip")

# Json library
include_directories(SYSTEM "libjsoncpp")

# Seqan library
include_directories(SYSTEM "libseqan")

# Seqan library
include_directories(SYSTEM "libnanoflann")

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

add_subdirectory(libdip)
add_subdirectory(libtsne)
add_subdirectory(libjsoncpp)
add_subdirectory(src)
add_subdirectory(test)

# installation 
message(STATUS "CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})
install(PROGRAMS ${CMAKE_BINARY_DIR}/acdc DESTINATION bin)
install(DIRECTORY share/ DESTINATION share)
