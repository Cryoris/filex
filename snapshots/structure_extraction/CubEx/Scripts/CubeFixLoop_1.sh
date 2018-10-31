#..modify this if necessary
datapath="../dataloc"
pixtablepath="../scipost"
sourcemaskname="SourceMask.fits"

#..below are "standard" file names, modify if necessary
cubename="DATACUBE_FINAL_WCS_"
pixtablename="PIXTABLE_REDUCED_0001.fits"
outname="DATACUBE_FINAL_WCS_FIX_"

#NB: SliceEdge-masks produced in the previous step will not be overwritten in this step by default (-edgemaskpix 0)

#----- nothing to modify after this line ------

#..check that arguments are present
if [ $# -lt 2 ]; then
    echo "Please provide the initial and final number of the cubes to reduce! E.g., CubeFixLoop_1.sh 1 18"
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

    CMD="CubeFix -cube \"$datapath/$cubename$i.fits\" -pixtable \"$pixtablepath/$i/$pixtablename\" -sourcemask \"$sourcemaskname\" -out \"$outname$i.fits\" -edgemaskpix 0 >& CubeFixLOG_1.$i &"
    echo $CMD
    eval $CMD

done



