#! /bin/bash

##
# Description:
# Run single_fits.sh for all snapshots in EAGLE_snapshots/RefL0025N0752
##

# One might want to modify this list, if not all elements are supposed to be processed
list=$(ls -lh /net/astrogate/export/astrodata/EAGLE_snapshots/RefL0025N0752/ | sed 1d | rev | cut -b1-12 | rev)

for el in $list; do
  # Skip directory if it already exists
  if [ -d "../snapshot_$el" ]; then
    echo "Directory snapshot_$el already exists, skipping.."
    continue
  fi
  echo "Running $el.."
  sh single_fits.sh $el
done
