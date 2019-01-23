#!/bin/bash

dir='parse_VolpoKhIsotherm'
rm -f parse_PE.out

while read f; do
  if [ -d $dir/$f ]; then
    echo "$f: computing PE"
    calPE.py $f coal -datapath ./$dir >> parse_PE.out
  else
    echo "$f: does not exist... still computing / error / nonporous / not submitted"
  fi
done < 3dN.list
