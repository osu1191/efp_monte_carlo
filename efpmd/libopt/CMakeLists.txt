set(SRC_OPT opt.c)
set(SRC_F
        blas.f
        lbfgsb.f
        linpack.f
        timer.f
)

enable_language(Fortran)
add_library(opt ${SRC_OPT} ${SRC_F} opt.h)
