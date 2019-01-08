cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# Maps to Visual Studio solution file (Tutorial.sln)
# The solution will have all targets (exe, lib, dll) 
# as Visual Studio projects (.vcproj)
project (ThermalFIST)

# The version number.
set (ThermalFIST_VERSION_MAJOR 0)
set (ThermalFIST_VERSION_MINOR 8)
set (ThermalFIST_VERSION_DEVEL 6)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/ThermalFISTConfig.h.in"
  "${PROJECT_BINARY_DIR}/include/ThermalFISTConfig.h"
  )

# Turn on the ability to create folders to organize projects (.vcproj)
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebug RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if(CMAKE_BUILD_TYPE STREQUAL "" AND NOT CMAKE_CXX_FLAGS MATCHES "-O[123]")
   message(STATUS "WARNING! It seems you are compiling without optimization. Please set CMAKE_BUILD_TYPE=Release for optimization.")
endif(CMAKE_BUILD_TYPE STREQUAL "" AND NOT CMAKE_CXX_FLAGS MATCHES "-O[123]")

# Set compiler flags and options. 
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#else()
#        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
# Here it is setting the Visual Studio warning level to 4
if(MSVC)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
else(MSVC)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif(MSVC)

if(MSVC)
else(MSVC)
#set(CMAKE_BUILD_TYPE Release)
endif(MSVC)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
#list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
#message ( ${CMAKE_PREFIX_PATH} )
#list(APPEND CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)

#find_package(ROOT COMPONENTS Core MathCore Matrix Minuit2)

#if(NOT STANDALONE_MINUIT)
#   set(STANDALONE_MINUIT 0)
#endif(NOT STANDALONE_MINUIT)

option (STANDALONE_MINUIT 
        "Build Minuit2 within project instead of linking one from ROOT distribution" OFF) 

string(COMPARE EQUAL "$ENV{ROOTSYS}" "" result)
if(result)
  set(ROOT_FOUND 0)
  set(STANDALONE_MINUIT 1)
else(result)
  set(ROOT_FOUND 1)
  find_library(MINUIT2_LIB NAMES Minuit2 libMinuit2.lib PATHS "$ENV{ROOTSYS}/lib/root/" "$ENV{ROOTSYS}/lib/")
  if(NOT MINUIT2_LIB)
    set(STANDALONE_MINUIT 1)
  endif()
endif(result)

#set(USE_MINUIT 1)
add_definitions(-DUSE_MINUIT)

#message (${ROOT_FOUND})

if (ROOT_FOUND)
add_definitions(-DUSE_ROOT)
else (ROOT_FOUND)
message(WARNING "ROOT not found! Minuit2 will be built as a standalone package")
endif (ROOT_FOUND)

OPTION (USE_OpenMP "Use OpenMP" OFF)
if(USE_OpenMP)
	find_package(OpenMP)
	if (OPENMP_FOUND)
		#message(STATUS "Found OpenMP!")
		add_definitions(-DUSE_OPENMP)
		SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	else (OPENMP_FOUND)
		message(WARNING "OpenMP not found! No parallelization.")
	endif (OPENMP_FOUND)
endif(USE_OpenMP)

# Command to output information to the console
# Useful for displaying errors, warnings, and debugging
message ("cxx Flags: " ${CMAKE_CXX_FLAGS})

# Sub-directories where more CMakeLists.txt exist

# Include directories of libraries
set (HRGBase_INCLUDE            "${PROJECT_SOURCE_DIR}/include/HRGBase")
set (HRGEV_INCLUDE              "${PROJECT_SOURCE_DIR}/include/HRGEV")
set (HRGVDW_INCLUDE             "${PROJECT_SOURCE_DIR}/include/HRGVDW")
set (HRGEventGenerator_INCLUDE  "${PROJECT_SOURCE_DIR}/include/HRGEventGenerator")
set (HRGFit_INCLUDE             "${PROJECT_SOURCE_DIR}/include/HRGFit")

#SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Main library
add_subdirectory(src/libraries)

# Gui
add_subdirectory(src/gui/QtThermalFIST)

# Sample routines
add_subdirectory(src/routines)

# build a CPack driven installer package
#include (InstallRequiredSystemLibraries)
#set (CPACK_RESOURCE_FILE_LICENSE  
#     "${CMAKE_CURRENT_SOURCE_DIR}/License.txt")
#set (CPACK_PACKAGE_VERSION_MAJOR "${CuteHRG_VERSION_MAJOR}")
#set (CPACK_PACKAGE_VERSION_MINOR "${CuteHRG_VERSION_MINOR}")
#include (CPack)



# Testing with GoogleTest (optional, but recommended)
option (INCLUDE_TESTS 
        "Unit testing with GoogleTest" OFF) 

if(INCLUDE_TESTS)

	# Download and unpack googletest at configure time
	configure_file(CMakeLists.txt.in ${PROJECT_BINARY_DIR}/googletest-download/CMakeLists.txt)
	execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
	  RESULT_VARIABLE result
	  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/googletest-download )
	if(result)
	  message(FATAL_ERROR "CMake step for googletest failed: ${result}")
	endif()
	execute_process(COMMAND ${CMAKE_COMMAND} --build .
	  RESULT_VARIABLE result
	  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/googletest-download )
	if(result)
	  message(FATAL_ERROR "Build step for googletest failed: ${result}")
	endif()

	# Prevent overriding the parent project's compiler/linker
	# settings on Windows
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

	# Add googletest directly to our build. This defines
	# the gtest and gtest_main targets.
	add_subdirectory(${PROJECT_BINARY_DIR}/googletest-src
					 ${PROJECT_BINARY_DIR}/googletest-build
					 EXCLUDE_FROM_ALL)

	# The gtest/gtest_main targets carry header search path
	# dependencies automatically when using CMake 2.8.11 or
	# later. Otherwise we have to add them here ourselves.
	if (CMAKE_VERSION VERSION_LESS 2.8.11)
	  include_directories("${gtest_SOURCE_DIR}/include")
	endif()

	enable_testing()

	add_subdirectory(test)

endif()