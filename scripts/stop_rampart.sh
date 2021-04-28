#!/bin/sh



lsof_out=$(lsof -n -i :3000 | grep -E "^node")

echo "lsof_out: ${lsof_out}"




if [ -n "${lsof_out}" ]; then
    echo "Rampart is running and will be killed"
    
    killpid=$(echo $lsof_out | awk '{ print $2 }')
    echo "the pid to kill is ${killpid}"

    kill -2 $killpid
else
    echo "Nothing to kill"
fi