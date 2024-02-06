#!/bin/bash
if [ $# -eq 0 ] 
then
    dcbtype=Release
elif [ "$1" = "-d" ]; then
    dcbtype=Debug
else
    echo "Incorrect args given. Use -d for Debug. Release will be used as default."
    dcbtype=Release
fi
echo "hej"
echo "cd ${dcbtype}/src/main/"
echo "./mainout"
cd ${dcbtype}/src/main/
./mainout
