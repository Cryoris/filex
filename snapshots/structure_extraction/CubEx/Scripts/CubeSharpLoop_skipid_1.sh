cubename="DATACUBE_FINAL_WCS_FIX_"
outname="DATACUBE_FINAL_WCS_FIX_SHARP_"
sourcemaskname="SourceMask.fits"
hsncubename="CombCube_0.fits"
objid=$3

#..check that arguments are present
if [ $# -lt 2 ]; then
    echo "Please provide the initial and final number of the cubes to reduce! E.g., CubeSharpLoop_1.sh 1 36"
    exit
fi

#..loop over files
for n in $(seq $1 $2); do

    #..adjust file number adding initial 0s
    if [ $n -lt 10 ]; then
	i="000$n"
    elif [ $n -lt 100 ]; then 
	i="00$n"
    elif [ $n -lt 1000 ]; then
	i="0$n"
    else
	i="$n"
    fi

   CMD="CubeSharp -cube \"$cubename$i.fits\" -out \"$outname$i.fits\" -sourcemask \"$sourcemaskname\" -hsncube \"$hsncubename\" -lsig 5 -skipobjid \"$objid\" >& CubeSharpLOG_1.$i &"
# CMD="CubeSharp -cube \"$cubename$i.fits\" -out \"$outname$i.fits\" -sourcemask \"$sourcemaskname\" >& CubeSharpLOG_1.$i &"
   echo $CMD
    eval $CMD

done


