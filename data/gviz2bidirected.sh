#!/bin/bash
for fin in data/zoo/gviz/*.gv ;
do
fout=$(dirname  $fin)"/cleaned/_"$(basename  $fin)
echo $fout
awk '!seen[$0]++' $fin | \
	sed "1s/graph/digraph/" | \
	awk '{split($0,a,/;|--/); if(a[2]=="") \
	{print a[1]}else{print a[2],"->",a[1]; print a[1],"->",a[2]}}' > $fout
done

#find  $PWD/data/zoo/gviz/cleaned -type -f> data/input.txt