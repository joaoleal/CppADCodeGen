#! /bin/sh -e
if [ ! -d tmp ] ; then mkdir tmp ; fi
./testCodeGen 
rm -r tmp;
