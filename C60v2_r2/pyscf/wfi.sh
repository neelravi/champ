#!/bin/bash

#first step convert from pyscf to trexio (converter avaialble just spherical)

python /home/landinez/codes/trexio_tools/src/converters/pyscf_to_trexio.py -c RHF.dump -o RHF_sph.hdf5


# second step trexio spherical to cartesian 
trexio convert-to -t cartesian -o RHF_cart.hdf5 RHF_sph.hdf5

##to hack trex2champ
# create a false gamess input file
touch fake.out

#third step from trexio to champ (regarding local converter inside current directory)


python trex2champ.py --trex RHF_cart.hdf5 --gamess fake.out --motype RHF --backend hdf5 --lcao --geometry --basis --pseudo  --determinants 
