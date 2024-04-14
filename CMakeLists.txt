# Specify the minimum version for CMake
cmake_minimum_required(VERSION 3.10)

# Project's name
project(miraxtract)

include_directories(include)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_STANDARD 99)

file( GLOB_RECURSE LIBXTRACT_SOURCE "src/*.c" )

add_library( xtract STATIC ${LIBXTRACT_SOURCE} )
target_compile_definitions( xtract PUBLIC -DUSE_OOURA -D_USE_MATH_DEFINES)


# Add executable target with source files listed in SOURCE_FILES variable
set(SOURCE_FILES tests/xttest_main.cpp tests/xttest_scalar.cpp tests/xttest_util.cpp)
add_executable(miraxtract_tests ${SOURCE_FILES})
target_link_libraries(miraxtract_tests xtract)

add_executable(miraxtract_example examples/simpletest/simpletest.cpp examples/simpletest/WaveFile.cpp)
target_link_libraries(miraxtract_example xtract)

file( GLOB_RECURSE LI_SOURCE "li/*.cpp" )
add_executable(li_test ${LI_SOURCE} examples/simpletest/WaveFile.cpp)
target_include_directories(li_test PRIVATE include)
target_link_libraries(li_test xtract)

# nanobind

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

# Create CMake targets for Python components needed by nanobind
find_package(Python 3.8 COMPONENTS Interpreter Development REQUIRED)

# Determine the nanobind CMake include path and register it
execute_process(
        COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NB_DIR)

message(STATUS "NanoBind Cmake directory: " ${NB_DIR})
list(APPEND CMAKE_PREFIX_PATH "${NB_DIR}")

# Import nanobind through CMake's find_package mechanism
find_package(nanobind CONFIG REQUIRED)

# We are now ready to compile the actual extension module
nanobind_add_module(
        miraxtract
        python/miraxtract.cpp
)
target_link_libraries(miraxtract PRIVATE xtract)

# Install directive for scikit-build
install(TARGETS miraxtract LIBRARY DESTINATION .)