#!/bin/sh
for i in results stepsize det diff; do
   cp $i.data "$i-$1.data"
done
