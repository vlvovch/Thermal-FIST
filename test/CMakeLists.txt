# Now simply link against gtest or gtest_main as needed. Eg
add_executable(test_IdealGasFunctions test_IdealGasFunctions.cpp)
target_link_libraries(test_IdealGasFunctions ThermalFIST GTest::gtest_main)
set_property(TARGET test_IdealGasFunctions PROPERTY FOLDER tests)


include(GoogleTest)

# IdealGasFunctions tests
gtest_discover_tests(test_IdealGasFunctions) 


# CPC tests

# Run cpc1HRGTDep
add_test(NAME RunCPC1a
         COMMAND cpc1HRGTDep 0
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Compare output with reference
add_test(NAME CompareCPC1a
         COMMAND ${CMAKE_COMMAND} -E compare_files ${CMAKE_BINARY_DIR}/cpc1.Id-HRG.TDep.out ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/cpc1.Id-HRG.TDep.out)

# Run cpc1HRGTDep
add_test(NAME RunCPC1b
         COMMAND cpc1HRGTDep 1
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Compare output with reference
add_test(NAME CompareCPC1b
         COMMAND ${CMAKE_COMMAND} -E compare_files ${CMAKE_BINARY_DIR}/cpc1.EV-HRG.TDep.out ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/cpc1.EV-HRG.TDep.out)

# Run cpc1HRGTDep
add_test(NAME RunCPC1c
         COMMAND cpc1HRGTDep 2
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Compare output with reference
add_test(NAME CompareCPC1c
         COMMAND ${CMAKE_COMMAND} -E compare_files ${CMAKE_BINARY_DIR}/cpc1.QvdW-HRG.TDep.out ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/cpc1.QvdW-HRG.TDep.out)

# Run cpc2HRGTDep a
add_test(NAME RunCPC2a
        COMMAND cpc2chi2 0
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Compare output with reference
add_test(NAME CompareCPC2a
        COMMAND ${CMAKE_COMMAND} -E compare_files ${CMAKE_BINARY_DIR}/cpc2.Id-HRG.ALICE2_76.chi2.TDep.out ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/cpc2.Id-HRG.ALICE2_76.chi2.TDep.out)

# Run cpc2HRGTDep b
add_test(NAME RunCPC2b
        COMMAND cpc2chi2 1
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})


add_test(NAME CompareCPC2b
        COMMAND ${CMAKE_COMMAND} -E compare_files ${CMAKE_BINARY_DIR}/cpc2.EV-HRG-TwoComponent.ALICE2_76.chi2.TDep.out ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/cpc2.EV-HRG-TwoComponent.ALICE2_76.chi2.TDep.out)

# Run cpc2HRGTDep c
add_test(NAME RunCPC2c
        COMMAND cpc2chi2 2
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_test(NAME CompareCPC2c
        COMMAND ${CMAKE_COMMAND} -E compare_files ${CMAKE_BINARY_DIR}/cpc2.EV-HRG-BagModel.ALICE2_76.chi2.TDep.out ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/cpc2.EV-HRG-BagModel.ALICE2_76.chi2.TDep.out)


# Run cpc2HRGTDep d
add_test(NAME RunCPC2d
        COMMAND cpc2chi2 3
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_test(NAME CompareCPC2d
        COMMAND ${CMAKE_COMMAND} -E compare_files ${CMAKE_BINARY_DIR}/cpc2.QvdW-HRG.ALICE2_76.chi2.TDep.out ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/cpc2.QvdW-HRG.ALICE2_76.chi2.TDep.out)

## Now add cpc3a
#add_test(NAME RunCPC3a
#        COMMAND cpc3chi2NEQ 0
#        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#
#add_test(NAME CompareCPC3a
#        COMMAND ${CMAKE_COMMAND} -E compare_files ${CMAKE_BINARY_DIR}/cpc3.EQ.chi2.out ${PROJECT_SOURCE_DIR}/src/examples/cpc/SampleOutput/cpc3.EQ.chi2.out)
#
## Now add cpc3b
#add_test(NAME RunCPC3b
#        COMMAND cpc3chi2NEQ 1
#        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
#
#add_test(NAME CompareCPC3b
#        COMMAND ${CMAKE_COMMAND} -E compare_files ${CMAKE_BINARY_DIR}/cpc3.NEQ.chi2.out ${PROJECT_SOURCE_DIR}/src/examples/cpc/SampleOutput/cpc3.NEQ.chi2.out)

# Now add cpc4a
add_test(NAME RunCPC4a
        COMMAND cpc4mcHRG 0
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

add_test(NAME CompareCPC4a
        COMMAND ${CMAKE_COMMAND} -E compare_files ${CMAKE_BINARY_DIR}/cpc4.analyt.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/cpc4.analyt.dat)