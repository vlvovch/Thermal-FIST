# Find the Qt libraries
find_package(Qt5Widgets)
find_package(Qt5PrintSupport)


if (Qt5Widgets_FOUND AND Qt5PrintSupport_FOUND)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

#include_directories ("${PROJECT_SOURCE_DIR}/HRGlib/include;${PROJECT_SOURCE_DIR}/EVHRGlib/include;${PROJECT_SOURCE_DIR}/HRGEventGeneratorlib/include;${PROJECT_SOURCE_DIR}/ThermalModelFitlib/include")
#include_directories ("${PROJECT_SOURCE_DIR}/include" "${PROJECT_BINARY_DIR}/include")

file (GLOB QTHRG_SOURCES 
      "*.h" , 
      "*.cxx" , 
	  "*.cpp")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS} ${Qt5PrintSupport_EXECUTABLE_COMPILE_FLAGS}")
	  
# Tell CMake to create the executable
add_executable(QtThermalFIST ${QTHRG_SOURCES})

if(WIN32)
  if(${ThermalFIST_VERSION_DEVEL} STREQUAL "0")
    set_property(TARGET QtThermalFIST PROPERTY OUTPUT_NAME "QtHRG-${ThermalFIST_VERSION_MAJOR}.${ThermalFIST_VERSION_MINOR}")
  else()
    set_property(TARGET QtThermalFIST PROPERTY OUTPUT_NAME "QtHRG-${ThermalFIST_VERSION_MAJOR}.${ThermalFIST_VERSION_MINOR}.${ThermalFIST_VERSION_DEVEL}")
  endif()
endif()

# Use the Widgets module from Qt 5.
target_link_libraries(QtThermalFIST ThermalFIST Qt5::Core Qt5::Widgets Qt5::PrintSupport)
set_property(TARGET QtThermalFIST PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")


set_property(TARGET QtThermalFIST PROPERTY FOLDER "GUI")

endif (Qt5Widgets_FOUND AND Qt5PrintSupport_FOUND)
