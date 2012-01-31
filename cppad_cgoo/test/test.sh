#! /bin/sh -e
if [ ! -d tmp ] ; then mkdir tmp ; fi
./testCGOO 
rm -r tmp;
