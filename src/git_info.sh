#!/bin/sh

#if [ -e 'git_info.inc' ]; then
#    rm git_info.inc
#fi

echo '! This include file is generated with git_info.sh' > git_info.inc

git log | head -3 | while read line
do
echo "write(unit,*) '$line'" >> git_info.inc
done

