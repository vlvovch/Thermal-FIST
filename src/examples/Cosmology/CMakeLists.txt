include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")


add_executable (example-CosmicTrajectory CosmicTrajectory.cpp)
target_link_libraries (example-CosmicTrajectory ThermalFIST)
set_property(TARGET example-CosmicTrajectory PROPERTY FOLDER "examples/Cosmology")
install (TARGETS example-CosmicTrajectory RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin/examples/Cosmology)

add_executable (TestEMMDerivatives TestEMMDerivatives.cpp)
target_link_libraries (TestEMMDerivatives ThermalFIST)
set_property(TARGET TestEMMDerivatives PROPERTY FOLDER "examples/Cosmology")
