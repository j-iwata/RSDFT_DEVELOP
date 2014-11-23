#!/bin/sh

START=`date +%s`
make -f Makefile.uspp
END=`date +%s`

echo `expr $END - $START`
