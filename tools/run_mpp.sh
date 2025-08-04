#!/bin/bash
f=$(find -name '*.yaml' | cut -d '/' -f 2)

N=$(sed '8q;d' $f | cut -d : -f 2 | cut -d ' ' -f 2)
x=$(sed '23q;d' $f | cut -d : -f 2 | cut -d [ -f 2 | cut -d ] -f 1)
s=$(echo "scale=10; $x" | bc)
a=$(echo "scale=10; $N * $s" | bc)
if [[ $(echo "$a < -1" | bc) == 1 ]]
then
  b=${a%.*}
  Ns=$(($b * -2))
  if [[ $Ns -lt 20 ]]
  then
    momentspp param=opt.bpp F=$f O=60 > /dev/null
  elif [[ $Ns -lt 70 ]]
  then 
    momentspp param=opt.bpp F=$f O=150 > /dev/null
  elif [[ $Ns -lt 150 ]]
  then    
    momentspp param=opt.bpp F=$f O=300 > /dev/null
  else 
    momentspp param=opt.bpp F=$f O=400 > /dev/null
  fi
else
  momentspp param=opt.bpp F=$f O=10 > /dev/null
fi

