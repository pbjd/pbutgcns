#!/bin/bash

trap "rm -rf $tmp; popd" EXIT SIGINT

# generate files in local temp
pushd $tmp
for u in $(awk '{print $1}' $utg)
do
    tigStore-adapter.py $cap $u || exit $?
    pbutgcns unitig_${u}.lay >> cns.fa || exit $?
done

mv cns.fa $cns || exit $?
