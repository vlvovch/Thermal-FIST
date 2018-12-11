set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY  "${CMAKE_BINARY_DIR}/lib")

# Base part of the library  
set(SRCS_HRGBase
HRGBase/IdealGasFunctions.cpp
HRGBase/NumericalIntegration.cpp
HRGBase/ThermalModelIdeal.cpp
HRGBase/ThermalModelBase.cpp
HRGBase/ThermalModelCanonical.cpp
HRGBase/ThermalModelCanonicalCharm.cpp
HRGBase/ThermalModelCanonicalStrangeness.cpp
HRGBase/ThermalParticle.cpp
HRGBase/ThermalParticleSystem.cpp
HRGBase/xMath.cpp
)

source_group("HRGBase\\Source Files" FILES ${SRCS_HRGBase})

set(HEADERS_HRGBase
${PROJECT_SOURCE_DIR}/include/HRGBase/IdealGasFunctions.h
${PROJECT_SOURCE_DIR}/include/HRGBase/BilinearSplineFunction.h
${PROJECT_SOURCE_DIR}/include/HRGBase/NumericalIntegration.h
${PROJECT_SOURCE_DIR}/include/HRGBase/SplineFunction.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalModelIdeal.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalModelBase.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalModelCanonical.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalModelCanonicalCharm.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalModelCanonicalStrangeness.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalParticle.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalParticleSystem.h
${PROJECT_SOURCE_DIR}/include/HRGBase/xMath.h
)	  

source_group("HRGBase\\Header Files" FILES ${HEADERS_HRGBase})

# Excluded-volume part of the library
set(SRCS_HRGEV
HRGEV/ThermalModelEVCrossterms.cpp
HRGEV/ThermalModelEVDiagonal.cpp
HRGEV/ExcludedVolumeHelper.cpp
HRGEV/ThermalModelEVCanonicalStrangeness.cpp
)

source_group("HRGEV\\Source Files" FILES ${SRCS_HRGEV})

set(HEADERS_HRGEV
${PROJECT_SOURCE_DIR}/include/HRGEV/ThermalModelEVCrossterms.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ThermalModelEVDiagonal.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ExcludedVolumeModel.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ExcludedVolumeHelper.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ThermalModelEVCanonicalStrangeness.h
)	

source_group("HRGEV\\Header Files" FILES ${HEADERS_HRGEV})

# van der Waals part of the library
set(SRCS_HRGVDW
HRGVDW/ThermalModelVDWFull.cpp
HRGVDW/ThermalModelVDWCanonicalStrangeness.cpp
)

source_group("HRGVDW\\Source Files" FILES ${SRCS_HRGVDW})

set(HEADERS_HRGVDW
${PROJECT_SOURCE_DIR}/include/HRGVDW/ThermalModelVDWFull.h
${PROJECT_SOURCE_DIR}/include/HRGVDW/ThermalModelVDWCanonicalStrangeness.h
)

source_group("HRGVDW\\Header Files" FILES ${HEADERS_HRGVDW})

# Fitter
file (GLOB HRGFIT_SOURCES_MINUIT
      "${PROJECT_SOURCE_DIR}/thirdparty/Minuit/src/*.cpp"
      "${PROJECT_SOURCE_DIR}/thirdparty/Minuit/src/*.cxx")

file (GLOB HRGFIT_HEADERS_MINUIT
      "${PROJECT_SOURCE_DIR}/thirdparty/Minuit/include/Minuit2/*.h")

file (GLOB HRGFIT_HEADERS_MINUIT_MATH
      "${PROJECT_SOURCE_DIR}/thirdparty/Minuit/include/Math/*.h")
	  
set(SRCS_HRGFit
HRGFit/ThermalModelFit.cpp
)

source_group("HRGFit\\Source Files" FILES ${SRCS_HRGFit})

set(HEADERS_HRGFit
${PROJECT_SOURCE_DIR}/include/HRGFit/ThermalModelFit.h
)

if (STANDALONE_MINUIT)
  set(SRCS_HRGFit
      ${SRCS_HRGFit}
      ${HRGFIT_SOURCES_MINUIT}
  )
  source_group("HRGFit\\Minuit2" FILES ${HRGFIT_SOURCES_MINUIT})
endif()


source_group("HRGFit\\Header Files" FILES ${HEADERS_HRGFit})	  

# Event generator part
set(SRCS_HRGEventGenerator
HRGEventGenerator/Acceptance.cpp
HRGEventGenerator/EventGeneratorBase.cpp
HRGEventGenerator/MomentumDistribution.cpp
HRGEventGenerator/ParticleDecays.cpp
HRGEventGenerator/RandomGenerators.cpp
HRGEventGenerator/SimpleEvent.cpp
HRGEventGenerator/SREventGenerator.cpp
HRGEventGenerator/SSHEventGenerator.cpp
)

source_group("HRGEventGenerator\\Source Files" FILES ${SRCS_HRGEventGenerator})

set(HEADERS_HRGEventGenerator
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/Acceptance.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/EventGeneratorBase.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/MomentumDistribution.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/ParticleDecays.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/RandomGenerators.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/SREventGenerator.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/SSHEventGenerator.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/SimpleEvent.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/SimpleParticle.h
)	

source_group("HRGEventGenerator\\Header Files" FILES ${HEADERS_HRGEventGenerator})

# Create named folders for the sources within the .vcproj
# Empty name lists them directly under the .vcproj
# source_group("HRG" FILES ${HRG_SOURCES})

# Properties->C/C++->General->Additional Include Directories
include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/thirdparty")
if (STANDALONE_MINUIT)
  include_directories("${PROJECT_SOURCE_DIR}/thirdparty/Minuit/include")
endif()

if (ROOT_FOUND AND NOT STANDALONE_MINUIT)
include_directories ("$ENV{ROOTSYS}/include/root") 
include_directories ("$ENV{ROOTSYS}/include") 
endif (ROOT_FOUND AND NOT STANDALONE_MINUIT)

# Set Properties->General->Configuration Type to Dynamic Library(.dll)
# Creates math.dll with the listed sources collected in the variables
# Also adds sources to the Solution Explorer
add_library(ThermalFIST ${SRCS_HRGBase} ${HEADERS_HRGBase} ${SRCS_HRGEV} ${HEADERS_HRGEV} ${SRCS_HRGVDW} ${HEADERS_HRGVDW} ${SRCS_HRGFit} ${HEADERS_HRGFit} ${SRCS_HRGEventGenerator} ${HEADERS_HRGEventGenerator})

if (ROOT_FOUND AND NOT STANDALONE_MINUIT)
find_library (MINUIT2_LIB NAMES Minuit2 libMinuit2.lib PATHS "$ENV{ROOTSYS}/lib/root/" "$ENV{ROOTSYS}/lib/")
message (${MINUIT2_LIB})
target_link_libraries(ThermalFIST ${MINUIT2_LIB})
endif (ROOT_FOUND AND NOT STANDALONE_MINUIT)

set_property(TARGET ThermalFIST PROPERTY FOLDER "libraries")

# Creates folder "libraries" and adds target project (math.vcproj)
# set_property(TARGET libHRG PROPERTY FOLDER "libraries")

#set_property(TARGET CuteHRG PROPERTY CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

target_include_directories(ThermalFIST PUBLIC "${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include" "${PROJECT_SOURCE_DIR}/thirdparty")

install (TARGETS ThermalFIST
         ARCHIVE DESTINATION ${PROJECT_BINARY_DIR}/lib)

		 
