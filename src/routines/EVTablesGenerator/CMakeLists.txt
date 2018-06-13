# Properties->C/C++->General->Additional Include Directories
include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")

set(SRCS
EVTablesGenerator.cpp
)

# Set Properties->General->Configuration Type to Application(.exe)
# Creates app.exe with the listed sources (main.cxx)
# Adds sources to the Solution Explorer
add_executable (EVTablesGenerator ${SRCS})

# Properties->Linker->Input->Additional Dependencies
target_link_libraries (EVTablesGenerator ThermalFIST)

# Creates a folder "executables" and adds target 
# project (app.vcproj) under it
set_property(TARGET EVTablesGenerator PROPERTY FOLDER "routines")
#set_property(TARGET EVTablesGenerator PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/routines")

# Adds logic to INSTALL.vcproj to copy app.exe to destination directory
install (TARGETS EVTablesGenerator
         RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin/routines)
		 