import sys
import uproot
import os
badfiles = list()

for files in os.scandir("/home/seehrar/projects/def-jmammei/seehrar/outputs/beam_develop-march2022"):
    file = uproot.open(files)
    if len(file.keys()) == 0:
        badfiles.append(files)
print(badfiles)
    