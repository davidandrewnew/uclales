#!/bin/bash

ln -sf ~/uclales/build/uclales uclales
ln -sf ~/uclales/cases/dcbl/NAMELIST NAMELIST
ln -sf ~/uclales/cases/dcbl/dcbl.pbs dcbl.pbs

rm dcbl.*.nc
rm *.err
rm *.out
rm *.rst
rm *.dcbl.*s

ls -la
