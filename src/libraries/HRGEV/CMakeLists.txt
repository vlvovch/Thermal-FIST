# Collect sources into the variable MATH_SOURCES without
# having to explicitly list each header and source file.
#
# CMake documentation states "We do not recommend using GLOB to collect a list of
# source files from your source tree. If no CMakeLists.txt file changes when a 
# source is added or removed then the generated build system cannot know when to 
# ask CMake to regenerate".
# Collect sources into the variable SIMPLE_FUNCTION_SOURCES

#file (GLOB EVHRG_SOURCES
#	  "*.cxx")
	  
#file (GLOB EVHRG_HEADERS
#      "${PROJECT_SOURCE_DIR}/include/HRGEV/*.h")
	  
set(SRCS
ThermalModelEVCrossterms.cxx
ThermalModelEVDiagonal.cxx
ExcludedVolumeHelper.cpp
ThermalModelEVCanonicalStrangeness.cxx
)

set(HEADERS
${PROJECT_SOURCE_DIR}/include/HRGEV/ThermalModelEVCrossterms.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ThermalModelEVDiagonal.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ExcludedVolumeModel.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ExcludedVolumeHelper.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ThermalModelEVCanonicalStrangeness.h
)	  

# Create named folders for the sources within the .vcproj
# Empty name lists them directly under the .vcproj
# source_group("HRG" FILES ${HRG_SOURCES})

# Properties->C/C++->General->Additional Include Directories
#include_directories ("include/;${PROJECT_SOURCE_DIR}/HRGlib/include")
include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")

if (ROOT_FOUND)
include_directories ("$ENV{ROOTSYS}/include/root") 
include_directories ("$ENV{ROOTSYS}/include") 
link_directories("$ENV{ROOTSYS}/lib/root/" "$ENV{ROOTSYS}/lib/")
endif (ROOT_FOUND)

# Set Properties->General->Configuration Type to Dynamic Library(.dll)
# Creates math.dll with the listed sources collected in the variables
# Also adds sources to the Solution Explorer
add_library(HRGEV ${SRCS} ${HEADERS})

if (ROOT_FOUND)
find_library (CORE_LIB NAMES Core libCore.lib PATHS "$ENV{ROOTSYS}/lib/root/" "$ENV{ROOTSYS}/lib/")
message (${CORE_LIB})
find_library (MATHCORE_LIB NAMES MathCore libMathCore.lib PATHS "$ENV{ROOTSYS}/lib/root/" "$ENV{ROOTSYS}/lib/")
message (${MATHCORE_LIB})
find_library (MATRIX_LIB NAMES Matrix libMatrix.lib PATHS "$ENV{ROOTSYS}/lib/root/" "$ENV{ROOTSYS}/lib/")
message (${MATRIX_LIB})
target_link_libraries(HRGEV ${CORE_LIB} ${MATHCORE_LIB} ${MATRIX_LIB} HRGBase)
else (ROOT_FOUND)
target_link_libraries(HRGEV HRGBase)
endif (ROOT_FOUND)

set_property(TARGET HRGEV PROPERTY FOLDER "libraries")

# Creates folder "libraries" and adds target project (math.vcproj)
# set_property(TARGET libHRG PROPERTY FOLDER "libraries")

install (TARGETS HRGEV
         ARCHIVE DESTINATION ${PROJECT_BINARY_DIR}/lib)
