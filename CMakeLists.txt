cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(acdc)

set(CMAKE_BUILD_TYPE Release)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -O2 -std=c++11")

# Boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost 1.55.0 COMPONENTS program_options system filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# yaml-cpp
add_subdirectory(libyamlcpp)
include_directories(SYSTEM "libyamlcpp/include")

# t-SNE
include_directories(SYSTEM "libtsne")

# Dip statistic
include_directories(SYSTEM "libdip")

# Seqan
include_directories(SYSTEM "libseqan")

# Nanoflann
include_directories(SYSTEM "libnanoflann")

# Catch
include_directories(SYSTEM "libcatch")

# Fastcluster
include_directories(SYSTEM "libfastcluster")

# Eigen 3
include_directories(SYSTEM "libeigen3")

# OpenMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

add_subdirectory(libdip)
add_subdirectory(libtsne)
add_subdirectory(libfastcluster)
add_subdirectory(src)
add_subdirectory(test)

# installation 
message(STATUS "CMAKE_INSTALL_PREFIX: " ${CMAKE_INSTALL_PREFIX})
install(PROGRAMS ${CMAKE_BINARY_DIR}/acdc DESTINATION bin)
install(DIRECTORY DESTINATION share/acdc)
install(DIRECTORY share/acdc/ DESTINATION share/acdc)
install(FILES share/acdc/tools/acdc-make-fastas DESTINATION bin PERMISSIONS OWNER_EXECUTE OWNER_WRITE OWNER_READ GROUP_EXECUTE GROUP_READ)
install(FILES share/acdc/tools/acdc-filter-fasta-by-name.awk DESTINATION bin PERMISSIONS OWNER_WRITE OWNER_READ GROUP_READ)
