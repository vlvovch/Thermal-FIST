# Find the Qt libraries
find_package(Qt5Widgets)
find_package(Qt5PrintSupport)


if (Qt5Widgets_FOUND AND Qt5PrintSupport_FOUND)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTORCC ON)

file (GLOB QTHRG_SOURCES 
      "*.h" , 
      "*.cxx" , 
	    "*.cpp")

file (GLOB QCUSTOMPLOT_SOURCE
      "${PROJECT_SOURCE_DIR}/thirdparty/QCustomPlot/*.cpp")

file (GLOB QCUSTOMPLOT_HEADER
      "${PROJECT_SOURCE_DIR}/thirdparty/QCustomPlot/*.h")

source_group("QCustomPlot" FILES ${QCUSTOMPLOT_HEADER} ${QCUSTOMPLOT_SOURCE})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS} ${Qt5PrintSupport_EXECUTABLE_COMPILE_FLAGS}")
	  
# Tell CMake to create the executable
add_executable(QtThermalFIST ${QTHRG_SOURCES} ${QCUSTOMPLOT_HEADER} ${QCUSTOMPLOT_SOURCE} qttfres.qrc qttficon.rc)

if(WIN32)
  if(${ThermalFIST_VERSION_DEVEL} STREQUAL "0")
    add_custom_command(TARGET QtThermalFIST 
                       POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:QtThermalFIST> $<TARGET_FILE:QtThermalFIST>"/../QtThermalFIST-${ThermalFIST_VERSION_MAJOR}.${ThermalFIST_VERSION_MINOR}.exe")
  else()
    add_custom_command(TARGET QtThermalFIST 
                       POST_BUILD
                       COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:QtThermalFIST> $<TARGET_FILE:QtThermalFIST>"/../QtThermalFIST-${ThermalFIST_VERSION_MAJOR}.${ThermalFIST_VERSION_MINOR}.${ThermalFIST_VERSION_DEVEL}.exe")
  endif()
endif()

# Use the Widgets module from Qt 5.
target_link_libraries(QtThermalFIST ThermalFIST Qt5::Core Qt5::Widgets Qt5::PrintSupport)

target_include_directories(QtThermalFIST PUBLIC "${PROJECT_SOURCE_DIR}/thirdparty/QCustomPlot")

set_property(TARGET QtThermalFIST PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set_property(TARGET QtThermalFIST PROPERTY FOLDER "GUI")

endif (Qt5Widgets_FOUND AND Qt5PrintSupport_FOUND)
