#!/bin/bash
echo $1 $2 $3

PYTHON=~/Desktop/Python-3.5.1/install/bin/python3

if [ $1 = 'e' ]; then
  rm -r __pycache__; $PYTHON exp.py
elif [ $1  = 'init' ]; then
  if [ $2  = 'd' ]; then
    LD_LIBRARY_PATH=/cac/u01/mfa51/Desktop/libpng-1.6.23/install/lib:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH
    echo "LD_LIBRARY_PATH= " $LD_LIBRARY_PATH
  fi
else
  echo "Argument did not match!"
fi