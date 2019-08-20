#!/bin/bash

ln -sf ~/uclales/build/uclales uclales
ln -sf ~/uclales/cases/arm/NAMELIST NAMELIST
ln -sf ~/uclales/cases/arm/ls_flux_in ls_flux_in
ln -sf ~/uclales/cases/arm/arm.pbs arm.pbs

rm core
rm arm.*.nc
rm *.err
rm *.out
rm *.rst
rm *.arm.*s

ls -la
