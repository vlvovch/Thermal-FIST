# Check C++17 support needed for Qt6
if (NOT USE_QT5)
  include(CheckCXXCompilerFlag)
  CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)

  if (COMPILER_SUPPORTS_CXX17)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    find_package(Qt6 COMPONENTS Core Widgets PrintSupport)
  else()
    message(WARNING "C++17 not supported. Cannot use Qt6.")
    set(CMAKE_CXX_STANDARD 11)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
  endif()
endif()

if (Qt6_FOUND)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(UseQt6 TRUE)
else()
  message(WARNING "Qt6 not found or not used. Reverting to Qt5...")
  find_package(Qt5 COMPONENTS Core Widgets PrintSupport)
  if (Qt5_FOUND)
    set(UseQt5 TRUE)
  else()
    message(STATUS "Qt5 and Qt6 not found. QtThermalFIST will not be built.")
    set(UseQt5 FALSE)
  endif()
endif()
#
#find_package(Qt5Widgets)
#find_package(Qt5PrintSupport)


#if (Qt5Widgets_FOUND AND Qt5PrintSupport_FOUND)'
if(UseQt5 OR UseQt6)

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

if(UseQt6)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt6Widgets_EXECUTABLE_COMPILE_FLAGS} ${Qt6PrintSupport_EXECUTABLE_COMPILE_FLAGS}")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${Qt5Widgets_EXECUTABLE_COMPILE_FLAGS} ${Qt5PrintSupport_EXECUTABLE_COMPILE_FLAGS}")
endif(UseQt6)

add_compile_definitions(QT_DISABLE_DEPRECATED_UP_TO=0x060700)

# Tell CMake to create the executable
add_executable(QtThermalFIST ${QTHRG_SOURCES} ${QCUSTOMPLOT_HEADER} ${QCUSTOMPLOT_SOURCE} qttfres.qrc qttficon.rc trajectoriestab.cpp)

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
if(UseQt6)
  target_link_libraries(QtThermalFIST ThermalFIST Qt6::Core Qt6::Widgets Qt6::PrintSupport)
else()
  target_link_libraries(QtThermalFIST ThermalFIST Qt5::Core Qt5::Widgets Qt5::PrintSupport)
endif(UseQt6)

target_include_directories(QtThermalFIST PUBLIC "${PROJECT_SOURCE_DIR}/thirdparty/QCustomPlot")

set_property(TARGET QtThermalFIST PROPERTY RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

set_property(TARGET QtThermalFIST PROPERTY FOLDER "GUI")

endif(UseQt5 OR UseQt6)
#endif (Qt5Widgets_FOUND AND Qt5PrintSupport_FOUND)
