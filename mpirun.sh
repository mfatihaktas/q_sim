#!/bin/bash
echo $1 $2 $3

NUM_NODE=1
NODES=""
BEGIN_NODE_ID=3
END_NODE_ID=$(($BEGIN_NODE_ID + $NUM_NODE))
for i in `seq $BEGIN_NODE_ID $END_NODE_ID`; do
  STR=""
  [ "$i" -lt 10 ] && STR+="0"
  NODES+="dell$STR$i"
  [ "$i" -lt $END_NODE_ID ] && NODES+=","
done
echo "NODES= $NODES"

MPIRUN=/opt/openmpi-1.7ft_b3/bin/mpirun # /usr/lib64/openmpi/bin/mpirun
PKILL=/usr/bin/pkill
PYTHON=/usr/bin/python3 # ~/Desktop/Python-3.5.1/install/bin/python3
DIR=/cac/u01/mfa51/Desktop/q_sim

NUM_REP=100 # 100
NUM_APP=10 # $((4*8)) # $((5*8)) # $((5*16))

if [ $1  == 'r' ]; then
  # eval NODES=\$$NODES
  # NODE_LIST=(${NODES//,/ } )
  # for i in `seq 1 $NUM_NODE`; do
  #   NODE=${NODE_LIST[$(($(($i - 1)) % ${#NODE_LIST[@]} )) ] }
  #   echo "run app $i on $NODE"
  #   LOG_F="log/app_"$i".log"
  #   $MPIRUN -n $NUM_APP -host $NODE \
  #     $PYTHON $DIR/app.py < /dev/null 2>&1 | tee $LOG_F & # > $LOG_F 2>&1 < /dev/null & #
  # done
  for i in `seq 1 $NUM_REP`; do
    echo "rep $i:"
    for j in `seq 1 $NUM_APP`; do
      LOG_F="log/app_"$i"_"$j".log"
      $PYTHON app.py < /dev/null 2>&1 | tee $LOG_F & # > $LOG_F 2>&1 < /dev/null &
    done
    wait
  done
elif [ $1  == 'cat' ]; then
  cat log/* >> perm_log/nrep_"$NUM_REP"_napp_"$NUM_APP".log
  # rm log/*.log
elif [ $1  == 'k' ]; then
  # NODE_LIST=(${NODES//,/ } )
  # for NODE in "${NODE_LIST[@]}"; do
  #   echo "$NODE $PKILL -f app.py"
  #   $MPIRUN -npernode 1 -host $NODE $PKILL -f app.py
  # done
  $PKILL -f app
else
  echo "Argument did not match!"
fi