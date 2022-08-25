set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY  "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY  "${CMAKE_BINARY_DIR}/lib")

# Properties->C/C++->General->Additional Include Directories
include_directories ("${PROJECT_SOURCE_DIR}/include" 
"${PROJECT_BINARY_DIR}/include" 
"${PROJECT_SOURCE_DIR}/thirdparty/MersenneTwister" 
"${EIGEN3_INCLUDE_DIR}")

# Base part of the library  
set(SRCS_HRGBase
HRGBase/Broyden.cpp
HRGBase/IdealGasFunctions.cpp
HRGBase/NumericalIntegration.cpp
HRGBase/ParticleDecay.cpp
HRGBase/ThermalModelIdeal.cpp
HRGBase/ThermalModelBase.cpp
HRGBase/ThermalModelCanonical.cpp
HRGBase/ThermalModelCanonicalCharm.cpp
HRGBase/ThermalModelCanonicalStrangeness.cpp
HRGBase/ThermalParticle.cpp
HRGBase/ThermalParticleSystem.cpp
HRGBase/Utility.cpp
HRGBase/xMath.cpp
)

source_group("HRGBase\\Source Files" FILES ${SRCS_HRGBase})

set(HEADERS_HRGBase
${PROJECT_SOURCE_DIR}/include/HRGBase/Broyden.h
${PROJECT_SOURCE_DIR}/include/HRGBase/IdealGasFunctions.h
${PROJECT_SOURCE_DIR}/include/HRGBase/BilinearSplineFunction.h
${PROJECT_SOURCE_DIR}/include/HRGBase/NumericalIntegration.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ParticleDecay.h
${PROJECT_SOURCE_DIR}/include/HRGBase/SplineFunction.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalModelIdeal.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalModelBase.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalModelCanonical.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalModelCanonicalCharm.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalModelCanonicalStrangeness.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalParticle.h
${PROJECT_SOURCE_DIR}/include/HRGBase/ThermalParticleSystem.h
${PROJECT_SOURCE_DIR}/include/HRGBase/xMath.h
${PROJECT_SOURCE_DIR}/include/HRGBase/Utility.h
)	  

source_group("HRGBase\\Header Files" FILES ${HEADERS_HRGBase})

# Excluded-volume part of the library
set(SRCS_HRGEV
HRGEV/ThermalModelEVCrossterms.cpp
HRGEV/ThermalModelEVCrosstermsLegacy.cpp
HRGEV/ThermalModelEVDiagonal.cpp
HRGEV/ExcludedVolumeHelper.cpp
HRGEV/ThermalModelEVCanonicalStrangeness.cpp
)

source_group("HRGEV\\Source Files" FILES ${SRCS_HRGEV})

set(HEADERS_HRGEV
${PROJECT_SOURCE_DIR}/include/HRGEV/ThermalModelEVCrossterms.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ThermalModelEVCrosstermsLegacy.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ThermalModelEVDiagonal.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ExcludedVolumeModel.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ExcludedVolumeHelper.h
${PROJECT_SOURCE_DIR}/include/HRGEV/ThermalModelEVCanonicalStrangeness.h
)	

source_group("HRGEV\\Header Files" FILES ${HEADERS_HRGEV})

# van der Waals part of the library
set(SRCS_HRGVDW
HRGVDW/ThermalModelVDW.cpp
HRGVDW/ThermalModelVDWCanonicalStrangeness.cpp
)

source_group("HRGVDW\\Source Files" FILES ${SRCS_HRGVDW})

set(HEADERS_HRGVDW
${PROJECT_SOURCE_DIR}/include/HRGVDW/ThermalModelVDW.h
${PROJECT_SOURCE_DIR}/include/HRGVDW/ThermalModelVDWCanonicalStrangeness.h
)

source_group("HRGVDW\\Header Files" FILES ${HEADERS_HRGVDW})

# Fitter
	  
set(SRCS_HRGFit
HRGFit/ThermalModelFit.cpp
HRGFit/ThermalModelFitParameters.cpp
)

source_group("HRGFit\\Source Files" FILES ${SRCS_HRGFit})

set(HEADERS_HRGFit
${PROJECT_SOURCE_DIR}/include/HRGFit/ThermalModelFit.h
${PROJECT_SOURCE_DIR}/include/HRGFit/ThermalModelFitParameters.h
${PROJECT_SOURCE_DIR}/include/HRGFit/ThermalModelFitQuantities.h
)


source_group("HRGFit\\Header Files" FILES ${HEADERS_HRGFit})	  

# Event generator part
set(SRCS_HRGEventGenerator
HRGEventGenerator/Acceptance.cpp
HRGEventGenerator/EventGeneratorBase.cpp
HRGEventGenerator/FreezeoutModels.cpp
HRGEventGenerator/MomentumDistribution.cpp
HRGEventGenerator/ParticleDecaysMC.cpp
HRGEventGenerator/RandomGenerators.cpp
HRGEventGenerator/SimpleEvent.cpp
HRGEventGenerator/SphericalBlastWaveEventGenerator.cpp
HRGEventGenerator/CylindricalBlastWaveEventGenerator.cpp
HRGEventGenerator/CracowFreezeoutEventGenerator.cpp
HRGEventGenerator/EventWriter.cpp
HRGEventGenerator/HepMCEventWriter.cpp
HRGEventGenerator/HypersurfaceSampler.cpp
)

source_group("HRGEventGenerator\\Source Files" FILES ${SRCS_HRGEventGenerator})

set(HEADERS_HRGEventGenerator
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/Acceptance.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/EventGeneratorBase.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/FreezeoutModels.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/MomentumDistribution.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/ParticleDecaysMC.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/RandomGenerators.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/SphericalBlastWaveEventGenerator.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/CylindricalBlastWaveEventGenerator.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/CracowFreezeoutEventGenerator.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/SimpleEvent.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/SimpleParticle.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/EventWriter.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/HepMCEventWriter.h
${PROJECT_SOURCE_DIR}/include/HRGEventGenerator/HypersurfaceSampler.h
)	

source_group("HRGEventGenerator\\Header Files" FILES ${HEADERS_HRGEventGenerator})

# Partial chemical equilibrium part
set(SRCS_HRGPCE
HRGPCE/ThermalModelPCE.cpp
)

source_group("HRGPCE\\Source Files" FILES ${SRCS_HRGPCE})

set(HEADERS_HRGPCE
${PROJECT_SOURCE_DIR}/include/HRGPCE/ThermalModelPCE.h
)	

source_group("HRGPCE\\Header Files" FILES ${HEADERS_HRGPCE})

if (NOT STANDALONE_MINUIT)
include_directories ("$ENV{ROOTSYS}/include/root") 
include_directories ("$ENV{ROOTSYS}/include") 
endif (NOT STANDALONE_MINUIT)

# Set Properties->General->Configuration Type to Dynamic Library(.dll)
# Creates math.dll with the listed sources collected in the variables
# Also adds sources to the Solution Explorer
add_library(ThermalFIST 
${SRCS_HRGBase} ${HEADERS_HRGBase} 
${SRCS_HRGEV} ${HEADERS_HRGEV} 
${SRCS_HRGVDW} ${HEADERS_HRGVDW} 
${SRCS_HRGFit} ${HEADERS_HRGFit} 
${SRCS_HRGEventGenerator} ${HEADERS_HRGEventGenerator}
${SRCS_HRGPCE} ${HEADERS_HRGPCE}
)

if (NOT STANDALONE_MINUIT)
find_library (MINUIT2_LIB NAMES Minuit2 libMinuit2.lib PATHS "$ENV{ROOTSYS}/lib/root/" "$ENV{ROOTSYS}/lib/")
message (${MINUIT2_LIB})
target_link_libraries(ThermalFIST ${MINUIT2_LIB})
else (NOT STANDALONE_MINUIT)
target_link_libraries(ThermalFIST Minuit2)
endif (NOT STANDALONE_MINUIT)

set_property(TARGET ThermalFIST PROPERTY FOLDER "libraries")

target_include_directories(ThermalFIST PUBLIC 
"${PROJECT_SOURCE_DIR}/include" 
"${PROJECT_BINARY_DIR}/include"
"${PROJECT_SOURCE_DIR}/thirdparty/MersenneTwister")

#install (TARGETS ThermalFIST
#         ARCHIVE DESTINATION ${PROJECT_BINARY_DIR}/lib)

		 
