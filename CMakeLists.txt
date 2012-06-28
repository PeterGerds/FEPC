CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(FEPC)

option(HAS_FFTW3 "Use fftw3" 0)

if (HAS_FFTW3)
    set(CMAKE_C_FLAGS "${CCMAKE_C_FLAGS} -DHAS_FFTW3")
endif(HAS_FFTW3)

include_directories(include)

ADD_SUBDIRECTORY(source)

enable_testing()
ADD_SUBDIRECTORY(test)


