#!/bin/bash
echo $1 $2 $3

SIM_NODES=""
BEGIN_SIMNODE_ID=1
END_SIMNODE_ID=$(($BEGIN_SIMNODE_ID + 10))
for i in `seq $BEGIN_SIMNODE_ID $END_SIMNODE_ID`; do
  STR=""
  [ "$i" -lt 10 ] && STR+="0"
  SIM_NODES+="dell$STR$i"
  [ "$i" -lt $END_SIMNODE_ID ] && SIM_NODES+=","
done
echo "SIM_NODES= $SIM_NODES"

MATLAB=/cac/u01/mfa51/Desktop/matlab_2016/install_/bin/./matlab
PYTHON=~/Desktop/Python-3.5.1/install/bin/python3
MPIRUN=/cac/u01/mfa51/Desktop/openmpi-1.10.2/install/bin/mpirun
PKILL=/usr/bin/pkill

rm -r __pycache__
if [ $1 = 'e' ]; then
  # LOG_F="$1.log"
  # $PYTHON exp.py < /dev/null 2>&1 | tee $LOG_F
  # rm *.png; $PYTHON exp.py
  $PYTHON simplex_exp.py
elif [ $1 = 'x' ]; then
  $PYTHON mixed_exp.py
elif [ $1 = 'm' ]; then
  # $MATLAB -r "run exp.m; quit"
  # rm *.png; $PYTHON mds_exp.py # > mds_exp.log
  if [ -z "$2" ]; then
    $PYTHON mds_exp.py
  else
    $PYTHON mds_exp.py & # > mds_exp.log
  fi
elif [ $1 = 'mm' ]; then
  $PYTHON mds_models.py
elif [ $1 = 's' ]; then
  # $PYTHON simplex_models.py
  $PYTHON sching_pareto.py
elif [ $1 = 'c' ]; then
  $PYTHON codes_stability.py
elif [ $1 = 'd' ]; then
  # rm *.png; $PYTHON det_models.py
  $PYTHON deneme.py
elif [ $1 = 'dep' ]; then
  $PYTHON deprecated.py
elif [ $1 = 'a' ]; then
  # $PYTHON arepeat_sim_components.py
  # $PYTHON arepeat_models.py
  $PYTHON arepeat_exp.py
elif [ $1 = 'g' ]; then
  $PYTHON google_data.py
elif [ $1 = 'f' ]; then
  $PYTHON fairness_sim.py
elif [ $1 = 'me' ]; then
  rm *.log *.png
  NUM_Q_MIN=3
  NUM_Q_MAX=10
  NODE_LIST=(${SIM_NODES//,/ } )
  for NUM_Q in `seq $NUM_Q_MIN $NUM_Q_MAX`; do
    NODE=${NODE_LIST[$(($(($NUM_Q - 1)) % ${#NODE_LIST[@]} )) ] }
      echo "run sim NUM_Q= $NUM_Q on $NODE"
      LOG_F="sim_"$(($NUM_Q - 1))".log"
      $MPIRUN -x LD_LIBRARY_PATH -n 1 -host $NODE \
        $PYTHON $PWD/exp.py --num_q=$NUM_Q < /dev/null 2>&1 | tee $LOG_F & # > $LOG_F 2>&1 < /dev/null & #
    done
elif [ $1 = 'p' ]; then
  $PYTHON
elif [ $1  = 'init' ]; then
  if [ $2  = 'd' ]; then
    # LD_LIBRARY_PATH=/cac/u01/mfa51/Desktop/libgd-2.2.2/install/lib:$LD_LIBRARY_PATH
    LD_LIBRARY_PATH=/cac/u01/mfa51/Desktop/libpng-1.6.23/install/lib:$LD_LIBRARY_PATH
    export LD_LIBRARY_PATH
    echo "LD_LIBRARY_PATH= " $LD_LIBRARY_PATH
  fi
elif [ $1 = 'k' ]; then
  NODE_LIST=(${SIM_NODES//,/ } )
  for NODE in "${NODE_LIST[@]}"; do
    $MPIRUN -npernode 1 -host $NODE $PKILL -f exp
  done
else
  echo "Argument did not match!"
fi
