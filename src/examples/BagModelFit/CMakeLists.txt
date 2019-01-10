# Properties->C/C++->General->Additional Include Directories
include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")

set(SRCS
BagModelFit.cpp
)

# Set Properties->General->Configuration Type to Application(.exe)
# Creates app.exe with the listed sources (main.cxx)
# Adds sources to the Solution Explorer
add_executable (BagModelFit ${SRCS})

# Properties->Linker->Input->Additional Dependencies
target_link_libraries (BagModelFit ThermalFIST)

# Creates a folder "executables" and adds target 
# project (app.vcproj) under it
set_property(TARGET BagModelFit PROPERTY FOLDER "examples")

# Adds logic to INSTALL.vcproj to copy app.exe to destination directory
install (TARGETS BagModelFit
         RUNTIME DESTINATION ${PROJECT_BINARY_DIR}/bin/examples)
		 