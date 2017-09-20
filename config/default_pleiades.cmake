# Pleiades
set(CMAKE_Fortran_COMPILER "ifort")
set(Fortran_COMPILER_WRAPPER ifort)

set(USER_Fortran_FLAGS "-traceback -r8 -ftz -extend_source -lmpi -mkl -w -vec-report0 -opt-report0 -axCORE-AVX2 -xSSE4.2")
set(USER_Fortran_FLAGS_RELEASE "-O3 -no-prec-div -xHOST -fp-model source")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

set(NETCDF_INCLUDE_DIR "/home1/danew/local/include")
set(NETCDF_LIB_1       "/home1/danew/local/lib/libnetcdff.a")
set(NETCDF_LIB_2       "/home1/danew/local/lib/libnetcdf.a")
set(HDF5_LIB_1         "/home1/danew/local/lib/libhdf5_hl.a")
set(HDF5_LIB_2         "/home1/danew/local/lib/libhdf5.a")
set(SZIP_LIB           "")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)
