include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")


add_executable (example-ThermodynamicsBQS ThermodynamicsBQS.cpp)
target_link_libraries (example-ThermodynamicsBQS ThermalFIST)
set_property(TARGET example-ThermodynamicsBQS PROPERTY FOLDER "examples/EoS")
install (TARGETS example-ThermodynamicsBQS RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin/examples/EoS)
		 
add_executable (example-SusceptibilitiesBQS SusceptibilitiesBQS.cpp)
target_link_libraries (example-SusceptibilitiesBQS ThermalFIST)
set_property(TARGET example-SusceptibilitiesBQS PROPERTY FOLDER "examples/EoS")
install (TARGETS example-SusceptibilitiesBQS RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin/examples/EoS)
