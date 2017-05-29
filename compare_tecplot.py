#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file1')
parser.add_argument('file2')
args = parser.parse_args()

maxDiff = 0.0
with open(args.file1) as f1:
    with open(args.file2) as f2:
        lines1 = f1.readlines()[5:]
        lines2 = f2.readlines()[5:]
        assert(len(lines1) == len(lines2))
        for i in range(0,len(lines1)):
            vals1 = lines1[i].split()
            vals2 = lines2[i].split()
            assert(len(vals1) == len(vals2))
            for j in range(0,len(vals1)):
                maxDiff = max(maxDiff, abs(float(vals1[j]) - float(vals2[j])))
print maxDiff
