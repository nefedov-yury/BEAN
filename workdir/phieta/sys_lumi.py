#!/usr/bin/env python3
#
# systematic uncertainties for luminosity (pb-1)
# ------------------------------------------------------------------

import math

# J/Psi scan data 2012 (BAM-00157)
lum12 = [ \
         (14.919,0.029,0.158), \
         (15.060,0.029,0.158), \
         (17.393,0.193,0.0), \
         (4.769,0.017,0.052), \
         (15.558,0.030,0.162), \
         (14.910,0.030,0.157), \
         (2.143,0.011,0.023), \
         (1.816,0.010,0.019), \
         (2.135,0.011,0.023), \
         (2.069,0.011,0.024), \
         (2.203,0.011,0.024), \
         (0.756,0.007,0.008), \
         (1.612,0.010,0.018), \
         (2.106,0.011,0.022), \
         (1.720,0.010,0.019), \
         (1.264,0.009,0.013), \
        ]

# J/Psi scan data (15th — 18th April) 2018
lum18 = [ \
         (2.50166, 0.0152412, 0.0), \
         (2.96535, 0.0164732, 0.0), \
         (5.10588, 0.0220927, 0.0), \
         (3.07306, 0.0171561, 0.0), \
         (1.70679, 0.0125341, 0.0), \
         (4.78731, 0.0215091, 0.0), \
         (5.75102, 0.0238353, 0.0), \
         (5.8398,  0.0242264, 0.0), \
        ]

# R-scan data 2015
# http://docbes3.ihep.ac.cn/cgi-bin/DocDB/ShowDocument?docid=653
lumR = [ \
        (105.53 ,0.025,0.897), \
        ( 15.960,0.010,0.140), \
        ( 16.046,0.010,0.091), \
        ( 15.849,0.010,0.106), \
        ( 17.315,0.011,0.119), \
        (126.21 ,0.029,0.896), \
       ]

del lumR[:-1] # only 3080!

# ------------------------------------------------------------------
# print luminosity and total error
def print_totlum():
    for l in lum12+lum18+lumR:
        err = math.sqrt(l[1]**2+l[2]**2)
        print(f'{l[0]:.3f} \\pm {err:.3f}')
#  print_totlum()

# print relative error for luminosity
def print_relerr_lum():
    for l in lum12+lum18+lumR:
        err = math.sqrt(l[1]**2+l[2]**2)
        rel = err/l[0]
        print(f'{rel:.2e} : {100*rel:.2f}%')
#  print_relerr_lum()

# print relative error as C-header
def prt_one_vec(lum,suf):
    print(f'   vector<double> RelSysLum{suf} = {{')
    print('      ',end='')
    for i,l in enumerate(lum):
        err = math.sqrt(l[1]**2+l[2]**2)
        rel = err/l[0]
        print(f'{rel:.2e}',end='')
        if i+1 == len(lum):
            print('\n',end='')
        else:
            print(', ',end='')
            if (i+1)%5==0:
                print('\n      ',end='')
    print('   };')

def print_header_lum():
    print('   // relative systematic uncertainties for luminosity')
    prt_one_vec(lum12,"12")
    prt_one_vec(lum18,"18")
    prt_one_vec(lumR,"R")

print_header_lum()
