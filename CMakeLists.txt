cmake_minimum_required(VERSION 3.21.1 FATAL_ERROR)

# Maps to Visual Studio solution file (Tutorial.sln)
# The solution will have all targets (exe, lib, dll) 
# as Visual Studio projects (.vcproj)
project (ThermalFIST)

# The version number.
set (ThermalFIST_VERSION_MAJOR 1)
set (ThermalFIST_VERSION_MINOR 5)
set (ThermalFIST_VERSION_DEVEL 0)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/ThermalFISTConfig.h.in"
  "${PROJECT_BINARY_DIR}/include/ThermalFISTConfig.h"
  )

# set(CMAKE_CXX_STANDARD 11 CACHE STRING "The C++ standard to use")
# set(CMAKE_CXX_STANDARD_REQUIRED ON)
if (CMAKE_CXX_STANDARD)
	message(STATUS "Using C++${CMAKE_CXX_STANDARD}")
	SET(CMAKE_CXX_STANDARD_REQUIRED ON)
	SET(CXX_COMPILER_FLAGS_CHECKED TRUE)
endif()

if (NOT CXX_COMPILER_FLAGS_CHECKED)
	include(CheckCXXCompilerFlag)
	CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
	CHECK_CXX_COMPILER_FLAG("-std=c++14" COMPILER_SUPPORTS_CXX14)
	CHECK_CXX_COMPILER_FLAG("-std=c++1y" COMPILER_SUPPORTS_CXX1Y)
	if(COMPILER_SUPPORTS_CXX17)
		#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
		set(CMAKE_CXX_STANDARD 17)
	elseif(COMPILER_SUPPORTS_CXX14 OR COMPILER_SUPPORTS_CXX1Y)
		#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
		set(CMAKE_CXX_STANDARD 14)
	else()
		# message(WARNING "The compiler ${CMAKE_CXX_COMPILER} has no C++14 or C++17 support. Please use a different C++ compiler.")
		#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
		set(CMAKE_CXX_STANDARD 11 REQUIRED)
	# elseif(COMPILER_SUPPORTS_CXX1Y)
	# 	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1y")
	# 	set(CMAKE_CXX_STANDARD 14)
	endif()
	SET(CXX_COMPILER_FLAGS_CHECKED TRUE)
#    SET(CXX_COMPILER_FLAGS_CHECKED TRUE CACHE BOOL "Checked CXX compiler flags" FORCE)
endif()

# Turn on the ability to create folders to organize projects (.vcproj)
# It creates "CMakePredefinedTargets" folder by default and adds CMake
# defined projects like INSTALL.vcproj and ZERO_CHECK.vcproj
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebug RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

if (NOT EIGEN3_INCLUDE_DIR)
   set(EIGEN3_INCLUDE_DIR "${PROJECT_SOURCE_DIR}/thirdparty/eigen-3.4.0" CACHE STRING
      "Choose the location (include directory) of the Eigen3 library to use."
      FORCE)
endif(NOT EIGEN3_INCLUDE_DIR)

if(CMAKE_BUILD_TYPE STREQUAL "" AND NOT CMAKE_CXX_FLAGS MATCHES "-O[123]")
   message(STATUS "WARNING! It seems you are compiling without optimization. Please set CMAKE_BUILD_TYPE=Release for optimization.")
endif(CMAKE_BUILD_TYPE STREQUAL "" AND NOT CMAKE_CXX_FLAGS MATCHES "-O[123]")

# Set compiler flags and options. 
#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
#if(COMPILER_SUPPORTS_CXX11)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#elseif(COMPILER_SUPPORTS_CXX0X)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#else()
#	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
#endif()
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98 -pedantic")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

# Here it is setting the Visual Studio warning level to 4
if(MSVC)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
else(MSVC)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")
endif(MSVC)


option (STANDALONE_MINUIT 
        "Build Minuit2 within project instead of linking one from ROOT distribution" OFF) 

option (USE_QT5 
		"Use Qt5 for GUI instead of Qt6" OFF) 

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

# Check whether the size of int is at least 32-bit
INCLUDE (CheckTypeSize)
check_type_size("int" SIZEOFINT BUILTIN_TYPES_ONLY)
if (${SIZEOFINT} LESS "4")
message(FATAL_ERROR "Size of int of at least 4 bytes is required. Your compiler has " ${SIZEOFINT} " bytes. Please get \
a different compiler.")
endif()
message(STATUS "Size of int is " ${SIZEOFINT} " bytes")

add_definitions(-DUSE_MINUIT)

if(NOT STANDALONE_MINUIT)
add_definitions(-DUSE_ROOT)
message(STATUS "Minuit2 from ROOT installation will be used")
else(NOT STANDALONE_MINUIT)
message(STATUS "Minuit2 library will be built as a standalone package")
add_subdirectory(thirdparty/Minuit2)
endif(NOT STANDALONE_MINUIT)

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

# The main library
add_subdirectory(src/library)

# GUI
add_subdirectory(src/gui/QtThermalFIST)

# Routines
add_subdirectory(src/routines)

# Examples
add_subdirectory(src/examples)

# NuSTEAM example
add_subdirectory(src/nusteam)


# Testing with GoogleTest (optional, but recommended)
option (INCLUDE_TESTS 
        "Unit testing with GoogleTest" OFF) 

if(INCLUDE_TESTS)

	# Fetch GoogleTest
	include(FetchContent)
	FetchContent_Declare(
			googletest
			# URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
			URL https://github.com/google/googletest/archive/refs/tags/release-1.12.1.zip # v1.12.1 (last version to support C++11)
	)
	# For Windows: Prevent overriding the parent project's compiler/linker settings
	set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
	FetchContent_MakeAvailable(googletest)

	enable_testing()

	add_subdirectory(test)

endif()
