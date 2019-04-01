#!/usr/bin/env python3

import glob
import sys
import subprocess
import itertools as it

strats = ["--gold", "--highest", "--largest"]

files = glob.glob("*.mbx")
combinations = []
for r in range(1, len(strats)+1):
    for comb in it.combinations(strats, r):
        combinations.append(list(comb))


for comb in combinations:
    total_time = 0
    for f in files:
        print(" "*30, end="\r")
        print(f, "with", comb, end="\r")
        sys.stdout.flush()
        cmd = ["ibexopt-sip", "-q", f] + comb
        completed = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
        time = float(completed.stdout.split()[1])
        total_time += time
    print("Total time for {} = {}".format(comb, total_time))
