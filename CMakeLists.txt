CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

PROJECT(FEPC)

include_directories(include)

ADD_SUBDIRECTORY(source)

enable_testing()
ADD_SUBDIRECTORY(test)


