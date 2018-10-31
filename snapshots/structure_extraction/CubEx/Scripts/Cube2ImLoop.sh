
if [ -z "$1" ]
  then
    echo "Please provide the name of the file with the list of input cubes to process! E.g., Cube2ImLoop.sh cube.list"
    exit
fi

for thisfile in `cat $1`; do
CMD="Cube2Im $thisfile"
echo $CMD
eval $CMD
done

