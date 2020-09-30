#!/bin/bash

cd Seafile/Aktuell/Masterarbeit/kac-zwanzig-heatbath

python3 ./plots.py

rsync -av --delete -e ssh ledwon@dong48.physik.hu-berlin.de:/users/stud/ledwon/Seafile/Aktuell/Masterarbeit/kac-zwanzig-heatbath/ ~/Seafile/Aktuell/Masterarbeit/kac-zwanzig-heatbath/

