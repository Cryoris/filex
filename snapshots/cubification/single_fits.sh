#!/bin/bash

##
# Description:
# Script to generate FITS images from EAGLE simulation
# EAGLE hdf5 --P2C--> chombo hdf5 --PlotStArt--> density_cube.fits & temperature_cube.fits
##

##
# Step 1: Create all local directories
##

eagle_dir="/net/astrogate/export/astrodata/EAGLE_snapshots/RefL0025N0752"
snaps_dir="/net/astrogate/export/astrodata/jgacon/filex/snapshots"
this_dir="$(pwd)"

if [ $# -lt 1 ]; then
  echo "Usage:   sh single_fits.sh <target_id>"
  echo "Example: sh single_fits.sh 010_z003p984"
  echo "List of snaps:"
  ls -lh $eagle_dir
  exit 1
fi

# Snapshot to use
target="$1"
echo "Using snap_$target."

# Snapshot directories: EAGLE (source) and local (target)
eagle_snap_dir="$eagle_dir/snapshot_$target"
local_snap_dir="$snaps_dir/snapshot_$target"

# Create & access local folder
mkdir -p $local_snap_dir
cd $local_snap_dir

# Create chombo & fits directories
chombo_dir="chombo"
fits_dir="fits3D"
cubex_dir="cubex"
mkdir -p $chombo_dir $fits_dir $cubex_dir

##
# Step 2: Execute P2C
##

cd $chombo_dir
res=256
p2c_input="$eagle_snap_dir/snap_$target.%d.hdf5"
p2c_output="snap_${target}_${res}_${res}_8_box"

#!! tmp skip since I need to re-run PlotStArt
#$this_dir/P2C/bin/P2C-1.2 -inp $p2c_input -out "$p2c_output.chombo.hdf5" -base_grid $res -lmax 0 | tee -a "$p2c_output.log"

##
# Step 3: Execute PlotStArt
##

cd ../$fits_dir
ps_input="./../$chombo_dir/$p2c_output.chombo.hdf5"
ps_density_output="snap_${target}_${res}_${res}_8_box_density"
ps_temp_output="snap_${target}_${res}_${res}_8_box_temp"

$this_dir/PlotStArt/bin/PlotStART-0.1.0.a -inp $ps_input -out "$ps_density_output.fits" -var uniform-density -ftype fits3D -lmax 0 | tee -a "$ps_density_output.log"
$this_dir/PlotStArt/bin/PlotStART-0.1.0.a -inp $ps_input -out "$ps_temp_output.fits" -var uniform-temperature -ftype fits3D -lmax 0 | tee -a "$ps_temp_output.log"

##
# Done!
##

echo "Done with $1!"
