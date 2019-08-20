#!/bin/bash

ln -sf ~/uclales/build/uclales uclales
ln -sf ~/uclales/cases/cuxart/NAMELIST NAMELIST
ln -sf ~/uclales/cases/cuxart/cuxart.pbs cuxart.pbs

rm cuxart.*.nc
rm *.err
rm *.out
rm *.rst
rm *.cuxart.*s

ls -la
