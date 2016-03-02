#!/bin/sh
SOURCE=$1
THREAD_NUM=$2
BIN=$3
JROOT=$4
LOG_TIME1=`date +%H:%M:%S`
$JROOT/bin/jellyfish count  -m 14 -o $BIN/output -c 40 -s 4G -t $THREAD_NUM $SOURCE
LOG_TIME2=`date +%H:%M:%S`
echo "kmercounting time: "from $LOG_TIME1 to $LOG_TIME2
$JROOT/bin/jellyfish dump -c -t -o $BIN/out $BIN/output
LOG_TIME3=`date +%H:%M:%S`
echo "dump time (txt transfer): "from $LOG_TIME2 to $LOG_TIME3
rm $BIN/output
