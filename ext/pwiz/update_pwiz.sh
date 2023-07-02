find . -type f  -print0 | while read -d $'\0' file
do
  echo $file
  cp ../pwiz_new/$file $file
done
