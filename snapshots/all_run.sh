#! /bin/bash

##
# Description:
# Run single_fits.sh for all snapshots in EAGLE_snapshots/RefL0025N0752
##

# One might want to modify this list, if not all elements are supposed to be processed
list=$(ls -lh | grep "snapshot_0*" | sed 1d | rev | cut -b1-12 | rev)

for el in $list; do
  echo "Running $el.."
  sh single_run.sh $el
done
