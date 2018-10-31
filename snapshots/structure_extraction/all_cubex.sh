#! /bin/bash

##
# Description:
# Run single_cubex.sh for all snapshots in ../ with parameters $cubex_par.
##

# Get list of target snaps
list=$(ls -lh .. | grep "snapshot" | rev | cut -b1-12 | rev)

# Parameter file
cubex_par="cbx.par"

# Extract identifying parameters from parameter file.
# These will be used in the folder naming, where the data is stored.
halo_thres=$(sed -nr "/Halo\_Threshold\ =/ s/.*Halo\_Threshold\ =\ ([0-9]+).*/\1/p" $cubex_par)
sn_thres=$(sed -nr "/SN\_Threshold\ =/ s/^SN\_Threshold\ =\ ([0-9]+).*/\1/p" $cubex_par)
nvox=$(sed -nr "/MinNVox\ =/ s/^MinNVox\ =\ ([0-9]+).*/\1/p" $cubex_par)

cubex_sub_dir="f${sn_thres}_h${halo_thres}_v${nvox}"
echo "Data will be stored in subdirectory: $cubex_sub_dir"

for el in $list; do
  target_dir="../snapshot_$el/cubex/$cubex_sub_dir"
  if [ -d $target_dir ]; then
    echo "$target_dir already exists. Skipping."
    continue
  else
    echo "Good to go, $target_dir does not exist yet!"
  fi
  echo "Running $el.."
  sh single_cubex.sh $el
done
