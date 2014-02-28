#!/bin/bash

# always remove the tmp dir after exiting
trap "rm -rf $tmp; popd" EXIT SIGINT


# generate files in local temp
pushd $tmp
for u in $(awk '{print $1}' $utg)
do
    set -x
    tigStore-adapter.py $cap $u || exit $?
    pbutgcns -j ${nproc-1} unitig_${u}.lay >> cns.fa || exit $?
    set +x
done

# sometimes, the utg file will be empty :-/, sooo no consensus
set -x
if [ -e cns.fa ]
then
    mv cns.fa $cns
else
    touch $cns
fi
