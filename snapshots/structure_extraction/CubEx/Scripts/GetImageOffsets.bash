SEx=/usr/bin/sextractor  # <---- change location of SExtractor here
GetImageOffset=GetImageOffsets.x # <---- change location of GetImageOffsets.x here, if necessary
# these are to trim edges for MUSE cubes
MinEDGE="50"  
MaxEDGE="250"  
# max offset in pixels
MAXOFFSET="50"

if [ -z "$1" ]
  then
    echo "Please provide the name of the file with the list of input cubes to process! E.g., GetImageOffsets.bash inplist.dat"
    exit
fi


# produce a default.param for SExtractor
echo NUMBER > default.param
echo X_IMAGE >> default.param
echo Y_IMAGE >> default.param
echo FLUX_APER >> default.param

rm -f cat.list
i=0
for thisfile in `cat $1`
do
    i=$((i+1))
    # produce a white-light image with Cube2Im
    imfile=IM.$i.fits
    CMD="Cube2Im -cube $thisfile -out $imfile"
    eval $CMD<<EOF 
EOF

    # run SExtractor on the image
    catfile=cat.${imfile}.sex
    echo $catfile >> cat.list
    $SEx $imfile -CATALOG_TYPE ASCII -DETECT_THRESH 2 -DETECT_MINAREA 30 -FILTER N -CATALOG_NAME $catfile   # <--- change here the SExtractor parameters if needed

done

# remove default.param file
rm -f default.param

refcube=`head -n 1 $1`
 
CMD="$GetImageOffset -catlist cat.list -out offsets.dat -min_edge $MinEDGE -max_edge $MaxEDGE -maxoffset $MAXOFFSET -wcscube $refcube -wcsout OUTPUT_WCS_ "
echo $CMD
eval $CMD <<EOF
EOF

echo "output files= offsets.dat"

exit
# cleanup
i=0
for thisfile in `cat $1`
do
   i=$((i+1))
   imfile=IM.$i.fits
   catfile=cat.${imfile}.sex
   rm -f $imfile
   rm -f $catfile
done
rm -f cat.list


