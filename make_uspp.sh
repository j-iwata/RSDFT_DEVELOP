#!/bin/sh

echo $1
START=`date +%s`
make lda0 -f Makefile.uspp -j$1
make all -f Makefile.uspp -j$1
END=`date +%s`
echo `expr $END - $START`
