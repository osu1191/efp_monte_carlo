set(CMAKE_Fortran_COMPILER "/apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/bin/mpif90")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "10.2.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/apps/spack/bell/apps/gcc/10.2.0-gcc-4.8.5-mgcgzho/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/apps/spack/bell/apps/gcc/10.2.0-gcc-4.8.5-mgcgzho/bin/gcc-ranlib")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/include;/apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/lib;/apps/spack/bell/apps/gcc/10.2.0-gcc-4.8.5-mgcgzho/lib/gcc/x86_64-pc-linux-gnu/10.2.0/finclude;/apps/spack/bell/apps/qt/5.12.5-gcc-4.8.5-crhkujb/include;/apps/spack/bell/apps/libtiff/4.0.10-gcc-4.8.5-ay66rcp/include;/apps/spack/bell/apps/netlib-lapack/3.8.0-gcc-10.2.0-u6wadoq/include;/apps/spack/bell/apps/openblas/0.3.8-gcc-10.2.0-yierafl/include;/apps/spack/bell/apps/hdf5/1.10.6-gcc-10.2.0-fa6zg3j/include;/apps/spack/bell/apps/boost/1.68.0-gcc-10.2.0-illfiwk/include;/apps/spack/bell/apps/libszip/2.1.1-gcc-10.2.0-ep7tepl/include;/apps/spack/bell/apps/zlib/1.2.11-gcc-10.2.0-vfja5fy/include;/apps/spack/bell-20231031/apps/mpc/1.1.0-gcc-4.8.5-falbgx3/include;/apps/spack/bell-20231031/apps/mpfr/3.1.6-gcc-4.8.5-cslxql4/include;/apps/spack/bell-20231031/apps/gmp/6.2.1-gcc-4.8.5-dse55eg/include;/apps/spack/bell/apps/anaconda/2020.11-py38-gcc-4.8.5-nhzhrm2/include;/apps/spack/bell/apps/gcc/10.2.0-gcc-4.8.5-mgcgzho/lib/gcc/x86_64-pc-linux-gnu/10.2.0/include;/usr/local/include;/apps/spack/bell/apps/gcc/10.2.0-gcc-4.8.5-mgcgzho/include;/apps/spack/bell/apps/gcc/10.2.0-gcc-4.8.5-mgcgzho/lib/gcc/x86_64-pc-linux-gnu/10.2.0/include-fixed;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "mpi_usempif08;mpi_usempi_ignore_tkr;mpi_mpifh;mpi;gfortran;m;gcc_s;gcc;quadmath;m;gcc_s;gcc;pthread;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/apps/spack/bell/apps/hwloc/1.11.11-gcc-10.2.0-mwxfudv/lib;/apps/spack/bell/apps/zlib/1.2.11-gcc-10.2.0-vfja5fy/lib;/usr/lib64;/apps/spack/bell/apps/openmpi/3.1.6-gcc-10.2.0-kn4ct52/lib;/apps/spack/bell/apps/netlib-lapack/3.8.0-gcc-10.2.0-u6wadoq/lib64;/apps/spack/bell/apps/gcc/10.2.0-gcc-4.8.5-mgcgzho/lib64;/apps/spack/bell/apps/gcc/10.2.0-gcc-4.8.5-mgcgzho/lib/gcc/x86_64-pc-linux-gnu/10.2.0;/lib64;/apps/spack/bell/apps/qt/5.12.5-gcc-4.8.5-crhkujb/lib;/apps/spack/bell/apps/libtiff/4.0.10-gcc-4.8.5-ay66rcp/lib;/apps/spack/bell/apps/openblas/0.3.8-gcc-10.2.0-yierafl/lib;/apps/spack/bell/apps/hdf5/1.10.6-gcc-10.2.0-fa6zg3j/lib;/apps/spack/bell/apps/boost/1.68.0-gcc-10.2.0-illfiwk/lib;/apps/spack/bell/apps/libszip/2.1.1-gcc-10.2.0-ep7tepl/lib;/apps/spack/bell/apps/gcc/10.2.0-gcc-4.8.5-mgcgzho/lib;/apps/spack/bell-20231031/apps/mpc/1.1.0-gcc-4.8.5-falbgx3/lib;/apps/spack/bell-20231031/apps/mpfr/3.1.6-gcc-4.8.5-cslxql4/lib;/apps/spack/bell-20231031/apps/gmp/6.2.1-gcc-4.8.5-dse55eg/lib;/apps/spack/bell/apps/anaconda/2020.11-py38-gcc-4.8.5-nhzhrm2/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
