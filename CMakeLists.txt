cmake_minimum_required(VERSION 3.16)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/target)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/target)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/target)

project(Mitsuba
  VERSION 0.6.0
  DESCRIPTION
    "Personal Compilable Version of the Mitsuba Renderer"
  LANGUAGES
    CXX C)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
message(STATUS "Build Mitsuba in ${CMAKE_BUILD_TYPE} mode")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set(CMAKE_CXX_FLAGS "-msse2 -fPIC")

add_compile_definitions(
  GLEW_MX
  MTS_DEBUGs
  SINGLE_PRECISION
  SPECTRUM_SAMPLES=3
  MTS_SSE MTS_HAS_COHERENT_RT)

include_directories(include)
include_directories(ext/openpgl/openpgl/include)
include_directories(ext/openvdb/nanovdb)
include_directories(ext/openvdb/openvdb)
include_directories(ext/tbb/include)
include_directories(ext/matplotlib-cpp)
include_directories(ext/nanoflann/include)
include_directories(/usr/include/python3.10)

# Load the configure script
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/data/cmake)
if(${CMAKE_SYSTEM_NAME} MATCHES Linux)
  include(MitsubaConfigLinux)
else()
  message(FATAL_ERROR "${CMAKE_SYSTEM_NAME} System is currently not supported")
endif()

# Build dependencies
add_subdirectory(ext/openpgl)
#add_subdirectory(ext/openvdb)
add_subdirectory(ext/openvdb/nanovdb/nanovdb)
add_subdirectory(ext/tbb)

# Build the source code
add_subdirectory(src)

# Copy the data files needed at runtime
file(COPY data/schema DESTINATION ${CMAKE_BINARY_DIR}/data REGEX ".*\\.xsd$|.*\\.xsl$")
file(COPY data/microfacet DESTINATION ${CMAKE_BINARY_DIR}/data REGEX ".*\\.dat$")
file(COPY data/ior DESTINATION ${CMAKE_BINARY_DIR}/data REGEX ".*\\.spd$")
