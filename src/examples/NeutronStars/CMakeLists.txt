include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")

add_executable (example-NeutronStars-CSHRG NeutronStars-CSHRG.cpp)
target_link_libraries (example-NeutronStars-CSHRG ThermalFIST)
set_property(TARGET example-NeutronStars-CSHRG PROPERTY FOLDER "examples/NeutronStars")
install (TARGETS example-NeutronStars-CSHRG RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin/examples/NeutronStars)
		 