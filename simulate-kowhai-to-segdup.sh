#!/bin/bash

echo "Kowhai-to-Segdup" > k2s.txt

for i in {11..100}
do
	~/git/kowhai/Release/kowhai --sim -nH 10 -nP 5 -nR 1 -rB 1.5 -pC 0.75 -pJ 0.5 --host-sets-rate >> k2s.txt
	cat ~/git/kowhai/for-segdup-from-kowhai.txt | xargs ~/git/segdup/Release/segdup -n 10000 -o trace >> k2s.txt
	mv segdup-trace.csv trace$i.csv 
done

