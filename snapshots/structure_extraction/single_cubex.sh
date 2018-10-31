#!/bin/bash

##
# Description:
# Run CubEx on the temperature restricted density cube with the parameters 
# from $cubex_par.
##

snaps_dir="/net/astrogate/export/astrodata/jgacon/filex/snapshots"

if [ $# -lt 1 ]; then
  echo "Usage:   sh run_cubex.sh <target_id>"
  echo "Example: sh run_cubex.sh 010_z003p984"
  echo "List of snaps:"
  ls -lh ..
  exit 1
fi

target="$1"

##
# Step 1: Create local directories
##

cubex_par="cbx.par"
local_snap_dir="$snaps_dir/snapshot_$target"

echo "Local snapshot directory: $local_snap_dir"
cd $local_snap_dir

# Create cubex subdirectory, if it doens't already exist
mkdir -p cubex
cd cubex
echo "Now we should be in the local cubex directory. We're in: $(pwd)"

halo_thres=$(sed -nr "/Halo\_Threshold\ =/ s/.*Halo\_Threshold\ =\ ([0-9]+).*/\1/p" $cubex_par)
sn_thres=$(sed -nr "/SN\_Threshold\ =/ s/^SN\_Threshold\ =\ ([0-9]+).*/\1/p" $cubex_par)
nvox=$(sed -nr "/MinNVox\ =/ s/^MinNVox\ =\ ([0-9]+).*/\1/p" $cubex_par)
cubex_sub_dir="f${sn_thres}_h${halo_thres}_v${nvox}"

echo "Extracted parameters: filament $sn_thres halo $halo_thres nvox $nvox"
echo "Data will be stored in $(pwd)/$cubex_sub_dir"

mkdir -p $cubex_sub_dir
cd $cubex_sub_dir

# Input file for cubex (must be called snap.fits)
ln -s "temperature_restricted_density.fits" "snap.fits"

##
# Step 2: Create constant 1 fits image to use for the signal-to-noise ratio
# 1.fits is filled with ones
##

./CubeArit snap.fits "*" 0. 0.fits
./CubeArit 0.fits "+" 1. 1.fits
rm 0.fits

##
# Step 3: Run CubEx and give the child a name
##

./CubEx/bin/CubEx $cubex_par
mv snap.Objects_Id.fits "ids_$1.fits"
