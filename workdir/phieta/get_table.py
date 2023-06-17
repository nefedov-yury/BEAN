#!/usr/bin/env python3
# print a summary table with fitting results

import os
import fnmatch
import re

Dir = 'mkk_inter/'
Start = re.compile('[ ]*Interference:')

files = []
for file in os.listdir(Dir):
    if fnmatch.fnmatch(file, 'mkk*.txt'):
        files.append(file)

files.sort()

print('\n')
print(
    ' Data     | Ndat| Nsb|     Nfit       | diff  |     Nphi       |')
print('', 62*'-')
for i, file in enumerate(files):
    #  if i == 5:
        #  break
    #  print(file)
    with open(Dir+file) as f:
        for line in f:
            if Start.match(line):
                title = line[Start.match(line).end():-1]
                #  print(title)

                # next line: Data
                line = f.readline()
                ln = line.strip().split()
                #  print(ln)
                Ndat = float(ln[2])
                Nsb = float(ln[4])
                #  print(Ndat,Nsb)

                # next three lines: Fit
                line = f.readline()
                #  ln=line.strip().split()
                #  print(ln)
                #  chi2 = float(ln[2])
                #  ndf  = int(ln[4])
                #  print(chi2,ndf)

                line = f.readline()
                #  ln=line.strip().split()
                #  print(ln)
                #  sigma = float(ln[1])
                #  errsm = float(ln[3])
                #  F     = float(ln[6])
                #  errF  = float(ln[8])
                #  print(sigma,errsm,F,errF)

                line = f.readline()
                ln = line.strip().split()
                #  print(ln)
                Nphi = float(ln[1])
                errNphi = float(ln[3])
                Nfit = float(ln[5])
                errNfit = float(ln[7])
                #  print(Nphi,errNphi,Nfit,errNfit)

                print(f' {file[4:-4]:8s} |'
                      f' {Ndat:3.0f} | {Nsb:2.0f} |'
                      f' {Nfit:5.1f} +/- {errNfit:4.1f} |'
                      f' {Ndat-Nfit:5.1f} |'
                      f' {Nphi:5.1f} +/- {errNphi:4.1f} |')
                break

print('', 62*'-', '\n')
