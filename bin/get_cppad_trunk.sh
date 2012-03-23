#! /bin/bash -e
# 
src="$HOME/cppad/trunk/$1"
if [ "$1" == "" ] || [ ! -e "$src" ]
then
	echo "usage: get_cppad_trunk.sh <destination>"
	echo "where $HOME/cppad/trunk/<destination> is the file to copy to here"
	exit 1
fi
destination="$1"
echo "cp $src $destination"
cp $src $destination

