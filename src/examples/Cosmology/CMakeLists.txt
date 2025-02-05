# Properties->C/C++->General->Additional Include Directories
include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")


add_executable (Cosmic-Trajectory CosmicTrajectory.cpp)
target_link_libraries (Cosmic-Trajectory ThermalFIST)
set_property(TARGET Cosmic-Trajectory PROPERTY FOLDER "examples/Cosmology")
install (TARGETS Cosmic-Trajectory RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin/examples/Cosmology)
		 