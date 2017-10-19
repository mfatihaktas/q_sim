#!/bin/bash
echo $1 $2 $3

if [ $1 = 'i' ]; then
  srun --partition=main --nodes=1 --ntasks=1 --cpus-per-task=1 --mem=2000 --time=03:00:00 --export=ALL --pty bash -i
else
  echo "Argument did not match!"
fi
