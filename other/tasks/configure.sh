#!/bin/bash
if [ $# -eq 0 ] 
then
    dcbtype=Release
elif [ "$1" = "-d" ]; then
    dcbtype=Debug
elif [ "$1" = "-r" ]; then
    dcbtype=Release
else
    echo "Incorrect args given. Use -d for Debug. Release will be used as default."
    dcbtype=Release
fi

echo "cmake -DCMAKE_BUILD_TYPE=${dcbtype} -S . -B ${dcbtype}/"
cmake -DCMAKE_BUILD_TYPE=${dcbtype} -S . -B ${dcbtype}/
