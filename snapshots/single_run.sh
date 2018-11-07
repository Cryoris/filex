#! /bin/bash

if [ $# -eq 0 ]; then
  echo "Usage: sh single_run.sh <target_id>"
  echo "Example: sh single_run.sh 028_z000p000"
  exit
fi

echo "Temperature restriction.."
cd temperature_restriction 
python restrict.py $1
cd ..

echo "Structure extraction.."
cd structure_extraction
sh single_cubex.sh $1
cd ..

echo "Conversion to csv.."
cd snapshot_$1/cubex/f8_h50_v100
python ../../../fits_to_csv/dump.py ids_$1.fits ids_$1.csv
