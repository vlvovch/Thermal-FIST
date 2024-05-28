set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/nusteam")


add_executable (NuSTEAM-HRG-EoS NuSTEAM-HRG-EoS.cpp)
target_link_libraries (NuSTEAM-HRG-EoS ThermalFIST)
set_property(TARGET NuSTEAM-HRG-EoS PROPERTY FOLDER "examples/NuSTEAM")
