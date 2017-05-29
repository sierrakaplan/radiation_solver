#!/usr/bin/env python

import fileinput

num = None
i = 0
xi = []
eta = []
w = []
tgt = None
for line in fileinput.input():
    line = line.strip()
    if num is None:
        num = int(line)
        i = num
        continue
    if i == num:
        i = 0
        if len(xi) == 0:
            tgt = xi
        elif len(eta) == 0:
            tgt = eta
        elif len(w) == 0:
            tgt = w
        else:
            break
    tgt.append(float(line))
    i += 1

print num

for i in range(0,num):
    if xi[i] > 0 and eta[i] > 0:
        print xi[i]
for i in range(0,num):
    if xi[i] > 0 and eta[i] < 0:
        print xi[i]
for i in range(0,num):
    if xi[i] < 0 and eta[i] > 0:
        print xi[i]
for i in range(0,num):
    if xi[i] < 0 and eta[i] < 0:
        print xi[i]

for i in range(0,num):
    if xi[i] > 0 and eta[i] > 0:
        print eta[i]
for i in range(0,num):
    if xi[i] > 0 and eta[i] < 0:
        print eta[i]
for i in range(0,num):
    if xi[i] < 0 and eta[i] > 0:
        print eta[i]
for i in range(0,num):
    if xi[i] < 0 and eta[i] < 0:
        print eta[i]

for i in range(0,num):
    if xi[i] > 0 and eta[i] > 0:
        print w[i]
for i in range(0,num):
    if xi[i] > 0 and eta[i] < 0:
        print w[i]
for i in range(0,num):
    if xi[i] < 0 and eta[i] > 0:
        print w[i]
for i in range(0,num):
    if xi[i] < 0 and eta[i] < 0:
        print w[i]
