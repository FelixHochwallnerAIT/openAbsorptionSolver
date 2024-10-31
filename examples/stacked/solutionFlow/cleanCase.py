#!/usr/bin/env python3

import os

os.system('rm -r 1* 2* 3* 4* 5* 6* 7* 8* 9* > /dev/null 2>&1')
os.system('rm -r postProcessing/ > /dev/null 2>&1')
os.system('rm log.* > /dev/null 2>&1')
os.system('rm -r constant/polyMesh > /dev/null 2>&1')
os.system('rm -r pythonCase.p > /dev/null 2>&1')
os.system('rm -r 0/C* > /dev/null 2>&1')
os.system('rm -r 0/phi > /dev/null 2>&1')
os.system('rm -r temp* > /dev/null 2>&1')
