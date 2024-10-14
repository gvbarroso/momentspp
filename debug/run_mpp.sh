#!/bin/bash

for((i=1; i<=6; i++))
do
  cd model_$i
  f=model_$i.yaml
  N=$(sed '6q;d' $f | cut -d : -f 2 | cut -d ' ' -f 2)
  x=$(sed '18q;d' $f | cut -d : -f 2 | cut -d [ -f 2 | cut -d ] -f 1)
  s=$(echo "scale=10; $x" | bc)
  a=$(echo "scale=10; $N * $s" | bc)
  b=${a%.*}
  z=$(($b * -1))
  w2=$(($z + $z))
  w3=$(($z + $z + $z))
  w4=$(($z + $z + $z + $z))
  momentspp params=opt.bpp F=$f O=$w2 M=$w2 > /dev/null
  momentspp params=opt.bpp F=$f O=$w3 M=$w3 > /dev/null
  momentspp params=opt.bpp F=$f O=$w4 M=$w4 > /dev/null
  cd ..
done
