add_executable(fepc_easy_test fepc_easy_test.c)
target_link_libraries(fepc_easy_test fepc m)

add_executable(stest stest.c)
target_link_libraries(stest fepc m)

add_executable(discont_test discont_test.c)
target_link_libraries(discont_test fepc m)

add_executable(sparsetest sparsetest.c)
target_link_libraries(sparsetest fepc m)

add_executable(fepc_cp_test fepc_cp_test.c)
target_link_libraries(fepc_cp_test fepc m)

if (HAS_FFTW3)
    target_link_libraries(fepc_easy_test fftw3)
    target_link_libraries(discont_test fftw3)
    target_link_libraries(sparsetest fftw3)
    target_link_libraries(fepc_cp_test fftw3)
    target_link_libraries(stest fftw3)
endif(HAS_FFTW3)

add_test(fepc_easy_test fepc_easy_test)
add_test(stest stest)
add_test(discont_test discont_test)
add_test(sparsetest sparsetest)
add_test(fepc_cp_test fepc_cp_test)

