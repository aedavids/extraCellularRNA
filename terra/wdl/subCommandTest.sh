
#set -x

# test concurrency using sub shell
# we want to let salmon run and collect stats at same time

# grouping using () cause command to run in sub shell

echo BEGIN

#subShellRet =  (
#time (
sh -c 'echo begin subshell;\
    for i in {1..5};\
    do\
        echo child $i;\
        sleep 1;\
    done;\
    \
    echo end subshell;\
    exit 123' &

for i in {1..10};
do

    echo parent $i
    sleep 2
done


echo "subShellRet: xxx$subShellRet xxx"

echo END
