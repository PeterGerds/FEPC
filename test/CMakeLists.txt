include_directories (${${PROJECT_NAME}_SOURCE_DIR}/src) 

link_directories (${${PROJECT_NAME}_BINARY_DIR}/src)


add_executable (fepc_easy_test fepc_easy_test.c)
target_link_libraries (fepc_easy_test fepc fftw3 m) 

add_executable (stest stest.c)
target_link_libraries (stest fepc fftw3 m)

add_executable (discont_test discont_test.c)
target_link_libraries (discont_test fepc fftw3 m)


add_executable (sparsetest sparsetest.c)
target_link_libraries (sparsetest fepc fftw3 m)


add_test(fepc_easy_test fepc_easy_test)
add_test(stest stest)
add_test(discont_test discont_test)
add_test(sparsetest sparsetest)