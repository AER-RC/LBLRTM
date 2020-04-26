#!/bin/bash

# stage inputs. it is assumed that the user follows the file naming
# convention that LBL expects (TAPE5, TAPE3, EMISSIVITY, REFLECTIVITY)
# doing volume mounts of cross sections until i make it a submodule
export INDIR=/LBLRTM/LBLRTM_In
cp -r $INDIR/xs .
cp $INDIR/TAPE5 .
cp $INDIR/TAPE3 .
cp $INDIR/FSCDXS .
#for INFILE in `ls $INDIR/*/`; do cp $INFILE .; done

# run LNFL
echo "Running LBL"
./lblrtm

# move outputs to directory that is volume mounted
export OUTDIR=/LBLRTM/LBLRTM_Out
mv TAPE* $OUTDIR
mv OD* $OUTDIR

exit
