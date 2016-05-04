enable_testing()

include_directories(${CHECK_INCLUDE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/pcapack/build/include})


set(TEST_SRC ${TEST_SRC} means.f90 sgns.f90)


add_executable(pcapack_test ${TEST_SRC})
target_link_libraries(pcapack_test pcapack ${CHECK_LIBRARIES})

add_test(statpack_test ${EXECUTABLE_OUTPUT_PATH}/pcapack_test)
