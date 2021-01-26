# Thunder
set(CMAKE_Fortran_COMPILER "ifort")
set(Fortran_COMPILER_WRAPPER mpif90)

set(USER_Fortran_FLAGS "-traceback -r8 -ftz -extend_source")
set(USER_Fortran_FLAGS_RELEASE "-O3 -no-prec-div -xHOST -fp-model source")
set(USER_Fortran_FLAGS_DEBUG "-fpe0 -O0 -g -check all -check nopointers -check noarg_temp_created")

set(NETCDF_INCLUDE_DIR  "/usr/local/other/netcdf4/4.7.4/intel/19.1.3.304/include")
set(NETCDF_LIB_1        "/usr/local/other/netcdf4/4.7.4/intel/19.1.3.304/lib/libnetcdff.a")
set(NETCDF_LIB_2        "/usr/local/other/netcdf4/4.7.4/intel/19.1.3.304/lib/libnetcdf.a")
set(PNETCDF_INCLUDE_DIR "/usr/local/other/netcdf3/pnetcdf/1.12.2/impi/19.1.3.304/include")
set(PNETCDF_LIB         "/usr/local/other/netcdf3/pnetcdf/1.12.2/impi/19.1.3.304/lib/libpnetcdf.a")
set(HDF5_LIB_1          "/usr/local/other/hdf5/1.13.0/intel/19.1.3.304/lib/libhdf5_hl.a")
set(HDF5_LIB_2          "/usr/local/other/hdf5/1.13.0/intel/19.1.3.304/lib/libhdf5.a")
set(SZIP_LIB            "")
set(LIBS ${NETCDF_LIB_1} ${NETCDF_LIB_2} ${PNETCDF_LIB} ${HDF5_LIB_1} ${HDF5_LIB_2} ${SZIP_LIB} m z curl)

