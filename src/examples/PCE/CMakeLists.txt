include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")

add_executable (example-PCE-Saha-LHC PCE-Saha-LHC.cpp)
target_link_libraries (example-PCE-Saha-LHC ThermalFIST)
set_property(TARGET example-PCE-Saha-LHC PROPERTY FOLDER "examples/EoS")
install (TARGETS example-PCE-Saha-LHC RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin/examples/PCE)
		 