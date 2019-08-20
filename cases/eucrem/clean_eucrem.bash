#!/bin/bash

ln -sf ~/uclales/build/uclales uclales
ln -sf ~/cases/eucrem/NAMELIST NAMELIST
ln -sf ~/cases/eucrem/ls_flux_in ls_flux_in
ln -sf ~/cases/eucrem/eucrem.pbs eucrem.pbs

rm core
rm eucrem.*.nc
rm *.err
rm *.out
rm *.rst
rm *.eucrem.*s

ls -la
