# Now simply link against gtest or gtest_main as needed. Eg
add_executable(test_IdealGasFunctions test_IdealGasFunctions.cpp)
target_link_libraries(test_IdealGasFunctions ThermalFIST GTest::gtest_main)
set_property(TARGET test_IdealGasFunctions PROPERTY FOLDER tests)

add_executable(test_QvdWHRG test_QvdWHRG.cpp)
target_link_libraries(test_QvdWHRG ThermalFIST GTest::gtest_main)
set_property(TARGET test_QvdWHRG PROPERTY FOLDER tests)

add_executable(test_CompareOutputs test_CompareOutputs.cpp)
target_link_libraries(test_CompareOutputs)
set_property(TARGET test_CompareOutputs PROPERTY FOLDER tests)

include(GoogleTest)

# IdealGasFunctions tests
gtest_discover_tests(test_IdealGasFunctions)

# QvdW-HRG stress tests
gtest_discover_tests(test_QvdWHRG) 


# CPC tests

# Testing with GoogleTest (optional, but recommended)
option (INCLUDE_ALL_TESTS
        "Include tests which may be `non-deterministic` among compilers/hardware (such as cpc1)" OFF)

if(INCLUDE_ALL_TESTS)
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
else()
    # Run cpc1HRGTDep
    add_test(NAME RunCPC1a
             COMMAND cpc1HRGTDep 0
             WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

    # Compare output with reference using test_CompareOutputs.cpp
    add_test(NAME CompareCPC1a
             COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/cpc1.Id-HRG.TDep.out ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/cpc1.Id-HRG.TDep.out)


        # Run cpc1HRGTDep
        add_test(NAME RunCPC1b
                 COMMAND cpc1HRGTDep 1
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
        # Compare output with reference using test_CompareOutputs.cpp
        add_test(NAME CompareCPC1b
                 COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/cpc1.EV-HRG.TDep.out ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/cpc1.EV-HRG.TDep.out)

        # Run cpc1HRGTDep
        add_test(NAME RunCPC1c
                 COMMAND cpc1HRGTDep 2
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
        # Compare output with reference using test_CompareOutputs.cpp
        add_test(NAME CompareCPC1c
                 COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/cpc1.QvdW-HRG.TDep.out ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/cpc1.QvdW-HRG.TDep.out)
endif(INCLUDE_ALL_TESTS)

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


# Thermodynamics (Id-HRG)
add_test(NAME RunThermodynamicsIdHRG1
        COMMAND example-ThermodynamicsBQS 0. 0. 0 ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/BQS_paramrange_1.in Thermodynamics-Id-HRG-output-1.dat
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareThermodynamicsIdHRG1
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/Thermodynamics-Id-HRG-output-1.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/Thermodynamics-Id-HRG-output-1.dat)

add_test(NAME RunThermodynamicsIdHRG2
        COMMAND example-ThermodynamicsBQS 0. 0. 0 ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/BQS_paramrange_2.in Thermodynamics-Id-HRG-output-2.dat
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareThermodynamicsIdHRG2
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/Thermodynamics-Id-HRG-output-2.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/Thermodynamics-Id-HRG-output-2.dat)

# Thermodynamics (QvdW-HRG)
add_test(NAME RunThermodynamicsQvdWHRG1
        COMMAND example-ThermodynamicsBQS 0.329 3.42 0 ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/BQS_paramrange_1.in Thermodynamics-QvdW-HRG-output-1.dat
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareThermodynamicsQvdWHRG1
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/Thermodynamics-QvdW-HRG-output-1.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/Thermodynamics-QvdW-HRG-output-1.dat)

add_test(NAME RunThermodynamicsQvdWHRG2
        COMMAND example-ThermodynamicsBQS 0.329 3.42 0 ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/BQS_paramrange_2.in Thermodynamics-QvdW-HRG-output-2.dat
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareThermodynamicsQvdWHRG2
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/Thermodynamics-QvdW-HRG-output-2.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/Thermodynamics-QvdW-HRG-output-2.dat)

# Thermodynamics (RG-HRG)
add_test(NAME RunThermodynamicsRGHRG1
        COMMAND example-ThermodynamicsBQS 0.329 3.42 1 ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/BQS_paramrange_1.in Thermodynamics-RG-HRG-output-1.dat
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareThermodynamicsRGHRG1
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/Thermodynamics-RG-HRG-output-1.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/Thermodynamics-RG-HRG-output-1.dat)

add_test(NAME RunThermodynamicsRGHRG2
        COMMAND example-ThermodynamicsBQS 0.329 3.42 1 ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/BQS_paramrange_2.in Thermodynamics-RG-HRG-output-2.dat
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareThermodynamicsRGHRG2
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/Thermodynamics-RG-HRG-output-2.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/Thermodynamics-RG-HRG-output-2.dat)

# Susceptibilities (Id-HRG)
add_test(NAME RunSusceptibilitiesIdHRG1
        COMMAND example-SusceptibilitiesBQS 0. 0. 0 ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/BQS_paramrange_1.in 0 Susceptibilities-Id-HRG-output-1.dat
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareSusceptibilitiesIdHRG1
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/Susceptibilities-Id-HRG-output-1.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/Susceptibilities-Id-HRG-output-1.dat)
add_test(NAME RunSusceptibilitiesIdHRG2
        COMMAND example-SusceptibilitiesBQS 0. 0. 0 ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/BQS_paramrange_2.in 1 Susceptibilities-Id-HRG-output-2.dat
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareSusceptibilitiesIdHRG2
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/Susceptibilities-Id-HRG-output-2.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/Susceptibilities-Id-HRG-output-2.dat)

# Susceptibilities (QvdW-HRG)
add_test(NAME RunSusceptibilitiesQvdWHRG1
        COMMAND example-SusceptibilitiesBQS 0.329 3.42 0 ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/BQS_paramrange_1.in 0 Susceptibilities-QvdW-HRG-output-1.dat
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareSusceptibilitiesQvdWHRG1
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/Susceptibilities-QvdW-HRG-output-1.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/Susceptibilities-QvdW-HRG-output-1.dat)

add_test(NAME RunSusceptibilitiesQvdWHRG2
        COMMAND example-SusceptibilitiesBQS 0.329 3.42 0 ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/BQS_paramrange_2.in 1 Susceptibilities-QvdW-HRG-output-2.dat
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareSusceptibilitiesQvdWHRG2
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/Susceptibilities-QvdW-HRG-output-2.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/Susceptibilities-QvdW-HRG-output-2.dat)

# Susceptibilities (RG-HRG)
add_test(NAME RunSusceptibilitiesRGHRG1
        COMMAND example-SusceptibilitiesBQS 0.329 3.42 1 ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/BQS_paramrange_1.in 0 Susceptibilities-RG-HRG-output-1.dat
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareSusceptibilitiesRGHRG1
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/Susceptibilities-RG-HRG-output-1.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/Susceptibilities-RG-HRG-output-1.dat)

add_test(NAME RunSusceptibilitiesRGHRG2
        COMMAND example-SusceptibilitiesBQS 0.329 3.42 1 ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/BQS_paramrange_2.in 1 Susceptibilities-RG-HRG-output-2.dat
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareSusceptibilitiesRGHRG2
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/Susceptibilities-RG-HRG-output-2.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/Susceptibilities-RG-HRG-output-2.dat)
        
# Neutron star EoS
add_test(NAME RunNeutronStarEoS
        COMMAND example-NeutronStars-CSHRG
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_test(NAME CompareNeutronStarEoS
        COMMAND test_CompareOutputs ${CMAKE_BINARY_DIR}/NSMatter-CSHRG-leptons-output.dat ${PROJECT_SOURCE_DIR}/test/ReferenceOutput/NSMatter-CSHRG-leptons-output.dat)