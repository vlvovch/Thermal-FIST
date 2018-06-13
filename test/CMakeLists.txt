# Now simply link against gtest or gtest_main as needed. Eg
add_executable(test_IdealGasFunctions test_IdealGasFunctions.cpp)
target_link_libraries(test_IdealGasFunctions ThermalFIST gtest_main)
set_property(TARGET test_IdealGasFunctions PROPERTY FOLDER tests)
add_test(NAME IdealGasFunctions COMMAND test_IdealGasFunctions)