# Properties->C/C++->General->Additional Include Directories
include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")



add_executable (cpc1HRGTDep cpc1-HRG-TDep.cpp)
target_link_libraries (cpc1HRGTDep ThermalFIST)
set_property(TARGET cpc1HRGTDep PROPERTY FOLDER "examples/cpc")

add_executable (cpc2chi2 cpc2-chi2-vs-T.cpp)
target_link_libraries (cpc2chi2 ThermalFIST)
set_property(TARGET cpc2chi2 PROPERTY FOLDER "examples/cpc")

add_executable (cpc3chi2NEQ cpc3-chi2NEQ.cpp)
target_link_libraries (cpc3chi2NEQ ThermalFIST)
set_property(TARGET cpc3chi2NEQ PROPERTY FOLDER "examples/cpc")

add_executable (cpc4mcHRG cpc4-mcHRG.cpp)
target_link_libraries (cpc4mcHRG ThermalFIST)
set_property(TARGET cpc4mcHRG PROPERTY FOLDER "examples/cpc")
