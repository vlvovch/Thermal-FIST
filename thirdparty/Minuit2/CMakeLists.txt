set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY  "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY  "${CMAKE_BINARY_DIR}/lib")

# Minuit2 headers 
file (GLOB HEADERS_MINUIT
      "Minuit/include/Minuit2/*.h")
      
file (GLOB HEADERS_MINUIT_MATH
      "Minuit/include/Math/*.h")

# Minuit2 sources 
file (GLOB SOURCES_MINUIT
      "Minuit/src/*.cpp"
      "Minuit/src/*.cxx")

add_library(Minuit2 ${HEADERS_MINUIT} ${HEADERS_MINUIT_MATH} ${SOURCES_MINUIT})

target_include_directories(Minuit2 PUBLIC "Minuit/include")

set_property(TARGET Minuit2 PROPERTY FOLDER "libraries")

#install (TARGETS Minuit2
#         ARCHIVE DESTINATION ${PROJECT_BINARY_DIR}/lib)

		 
