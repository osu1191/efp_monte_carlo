set(SRC_FF ff.c)
#include_directories(../include)
add_library(ff STATIC ${SRC_FF} ff.h)
