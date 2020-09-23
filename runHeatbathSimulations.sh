#!/bin/sh
python3 kac-zwanzig.py

for i in 1 2 3 4 5
do
  nohup python3 ./kac-zwanzig.py &
done
