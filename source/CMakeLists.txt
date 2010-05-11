file(GLOB ALL_SOURCE_FILES *.c)

# Find all header files in the include directory
file(GLOB ALL_HEADER_FILES ${CMAKE_SOURCE_DIR}/include/*.h)

ADD_LIBRARY(fepc SHARED ${ALL_SOURCE_FILES} ${ALL_HEADER_FILES})
