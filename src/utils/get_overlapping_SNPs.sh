#!/bin/bash

FILE="./READ_output_ordered"
if [ $# -ne 0 ]; then
    FILE=$1
fi

tail -n+2 $FILE | awk '{a[$1]+=$10}END{for (i in a) print i, a [i]}' | sort | column -t
