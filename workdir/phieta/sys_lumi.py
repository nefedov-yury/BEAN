#!/usr/bin/env python3
#
# systematic uncertainties for luminosity (pb-1)
# ------------------------------------------------------------------

import math

# J/Psi scan data 2012 (BAM-00157)
#  Int. L (pb−1)
#  14.919±0.029±0.158
#  15.060±0.029±0.158
#  17.393±0.035±0.190
 #  4.769±0.017±0.052
#  15.558±0.030±0.162
#  14.910±0.030±0.157
 #  2.143±0.011±0.023
 #  1.816±0.010±0.019
 #  2.135±0.011±0.023
 #  2.069±0.011±0.024
 #  2.203±0.011±0.023
 #  0.756±0.007±0.008
 #  1.612±0.010±0.018
 #  2.106±0.011±0.022
 #  1.720±0.010±0.019
 #  1.264±0.009±0.013

lum12 = [ \
         (14.919, 0.029, 0.158), \
         (15.060, 0.029, 0.158), \
         (17.393, 0.033, 0.190), \
         ( 4.769, 0.017, 0.052), \
         (15.558, 0.030, 0.162), \
         (14.910, 0.030, 0.157), \
         ( 2.143, 0.011, 0.023), \
         ( 1.816, 0.010, 0.019), \
         ( 2.135, 0.011, 0.023), \
         ( 2.069, 0.011, 0.024), \
         ( 2.203, 0.011, 0.024), \
         ( 0.756, 0.007, 0.008), \
         ( 1.612, 0.010, 0.018), \
         ( 2.106, 0.011, 0.022), \
         ( 1.720, 0.010, 0.019), \
         ( 1.264, 0.009, 0.013), \
        ]

# J/Psi scan data 2018
  #  Lγγ, pb−1
#  2.467 ± 0.017
#  2.918 ± 0.019
#  4.978 ± 0.028
#  3.103 ± 0.020
#  1.681 ± 0.013
#  4.656 ± 0.027
#  5.638 ± 0.031
#  5.716 ± 0.031

lum18 = [ \
         (2.467, 0.017, 0.0), \
         (2.918, 0.019, 0.0), \
         (4.978, 0.028, 0.0), \
         (3.103, 0.020, 0.0), \
         (1.681, 0.013, 0.0), \
         (4.656, 0.027, 0.0), \
         (5.638, 0.031, 0.0), \
         (5.716, 0.031, 0.0), \
        ]

# R-scan data 2015
# http://docbes3.ihep.ac.cn/cgi-bin/DocDB/ShowDocument?docid=653
lumR = [ \
        (105.530, 0.025, 0.897), \
        ( 15.960, 0.010, 0.140), \
        ( 16.046, 0.010, 0.091), \
        ( 15.849, 0.010, 0.106), \
        ( 17.315, 0.011, 0.119), \
        (126.210, 0.029, 0.896), \
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
