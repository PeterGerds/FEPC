file(GLOB ALL_SOURCE_FILES *.c)

# Find all header files in the include directory
file(GLOB ALL_HEADER_FILES ${CMAKE_SOURCE_DIR}/include/*.h)

LIST(REMOVE_ITEM ALL_SOURCE_FILES ${CMAKE_CURRENT_SOURCE_DIR}/gauss.c)
LIST(REMOVE_ITEM ALL_HEADER_FILES ${CMAKE_SOURCE_DIR}/include/gauss.h)
	
ADD_LIBRARY(fepc SHARED ${ALL_SOURCE_FILES} ${ALL_HEADER_FILES})

if (EXISTS ${CMAKE_SOURCE_DIR}/include/gauss.h)
	ADD_LIBRARY(gauss SHARED gauss.c ${CMAKE_SOURCE_DIR}/include/gauss.h)
endif()

