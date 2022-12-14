#!/usr/bin/env python
#
# calculate total systematic uncertainties

from math import *

# ------------------------------------------------------------------
# Luminosity
# 3050-3120: (BAM-00157)
lum = [ \
 (14.919,0.029,0.158), \
 (15.060,0.029,0.158), \
 (30.942,0.044,0.338), \
 (4.769,0.017,0.052), \
 (15.558,0.030,0.162), \
 (14.910,0.030,0.157), \
 (2.143,0.011,0.023), \
 (1.816,0.010,0.019), \
 (2.135,0.011,0.023), \
 (2.069,0.011,0.024), \
 (2.203,0.011,0.023), \
 (0.756,0.007,0.008), \
 (1.612,0.010,0.018), \
 (2.106,0.011,0.022), \
 (1.720,0.010,0.019), \
 (1.264,0.009,0.013), \
      ]
# 3080new: http://docbes3.ihep.ac.cn/cgi-bin/DocDB/ShowDocument?docid=653
lum.append((126.21,0.029,0.896))

#  print lum

# print luminosity and error
def print_totlum():
    for l in lum:
        err=sqrt(l[1]**2+l[2]**2)
        print '{L:.3f} +/- {E:.3f}'.format(L=l[0],E=err)
#  print_totlum()

# print relative errors for luminosity
def print_relerr_lum():
    for l in lum:
        stat=100*l[1]/l[0]
        sys=100*l[2]/l[0]
        tot=100*sqrt(l[1]**2+l[2]**2)/l[0]
        print 'stat: {:.2f} sys: {:.2f} tot:{:.2f}'.format(stat,sys,tot)
#  print_relerr_lum()

err_lum=[100*sqrt(l[1]**2+l[2]**2)/l[0] for l in lum]
#  print 'err_lum=',err_lum
# ------------------------------------------------------------------
# Track reconstruction: 2 * 1%
err_tr=2.

# Photon reconstruction: 2 * 1%
err_ph=.2

# Branching fractions (PDG):
err_br=1.1

# ------------------------------------------------------------------
# Selection:
err_chi2=1.9
err_chi2E= [ 6.87, 2.57, 1.90, 18.77, 1.60,\
             1.51, 3.95, 3.07,  0.85, 1.02,\
             0.75, 1.45, 1.38,  2.09, 3.49,\
             8.14, 1.90 ]

err_mkk=2.4
err_mkkE= [ 1.83, 6.20,  1.83, 1.85, 4.76,\
            1.24, 4.33,  2.00, 0.73, 2.35,\
            2.78, 1.71, 10.41, 1.79, 8.51,\
            2.12, 0.85 ]

err_mgg=5.8
err_mggE= [ 4.71, 1.19,  6.84, 30.19, 0.97,\
            2.68, 6.45,  2.15,  1.71, 0.25,\
            2.80, 2.85, 14.44,  9.06, 8.07,\
           27.46, 5.80 ]

# ------------------------------------------------------------------
# MCGPJ:
err_mcgpj= [ 0.82, 1.70, 2.63, 3.23, 3.51,\
             2.91, 0.54, 0.15, 0.07, 0.06,\
             0.14, 0.30, 1.25, 1.97, 2.86,\
             3.08, 2.63 ]

# ------------------------------------------------------------------
# Sigma:
sig= [ 2.52517e+01, 2.80780e+01, 3.25385e+01, 1.68147e+01, 3.04337e+01,\
       4.08263e+01, 2.23957e+02, 7.54028e+02, 2.23143e+03, 3.58859e+03,\
       1.72718e+03, 4.52804e+02, 1.10674e+02, 6.47354e+01, 8.99491e+01,\
       4.05821e+01, 2.49830e+01 ]
errsig= [\
        5.26534e+00, 5.50655e+00, 5.50001e+00, 8.89751e+00, 5.91524e+00,\
        6.53745e+00, 3.95903e+01, 7.98528e+01, 1.28580e+02, 1.63469e+02,\
        1.09425e+02, 9.44162e+01, 3.68914e+01, 2.38559e+01, 2.84444e+01,\
        2.34301e+01, 1.81638e+00\
        ]

# ------------------------------------------------------------------
def print_totuns():
    str='  double err_sys[] = {'
    for i in range(17):
        err=sqrt(\
            err_lum[i]**2 +\
            err_tr**2 +\
            err_ph**2 +\
            err_br**2 +\
            err_chi2**2 +\
            err_mkk**2 +\
            err_mgg**2 +\
            err_mcgpj[i]**2\
                )
#          print '{E:.2f}%'.format(E=err)
        str = str + ' {E:.2f},'.format(E=err)
        if (i+1)%5==0: str = str + '\n' + ' '*22
    str = str[:-1] + ' };'
    print str

#  print_totuns()

# ------------------------------------------------------------------
def print_crosssection():
    for i in range(17):
        err=sqrt(\
            err_lum[i]**2 +\
            err_tr**2 +\
            err_ph**2 +\
            err_br**2 +\
            err_chi2**2 +\
            err_mkk**2 +\
            err_mgg**2 +\
            err_mcgpj[i]**2\
                )
        Sig=sig[i]*1e-3
        Stat=errsig[i]*1e-3
        Sys=Sig*err*1e-2
        print ' {0:.3f} $\\pm$ {1:.3f} $\\pm$ {2:.3f} \\\\'\
                .format(Sig,Stat,Sys)

print_crosssection()

# ------------------------------------------------------------------
def print_totuns_noE():
    str='  double err_sys[] = {'
    for i in range(17):
        err=sqrt(\
            err_lum[i]**2 +\
            err_mcgpj[i]**2\
                )
        str = str + ' {E:.2f},'.format(E=err)
        if (i+1)%5==0: str = str + '\n' + ' '*22
    str = str[:-1] + ' };'
    print str

#  print_totuns_noE()

# ------------------------------------------------------------------
def print_totuns_fullE():
    str='  double err_sys[] = {'
    for i in range(17):
        err=sqrt(\
            err_lum[i]**2 +\
            err_chi2E[i]**2 +\
            err_mkkE[i]**2 +\
            err_mggE[i]**2 +\
            err_mcgpj[i]**2\
                )
        str = str + ' {E:.2f},'.format(E=err)
        if (i+1)%5==0: str = str + '\n' + ' '*22
    str = str[:-1] + ' };'
    print str

#  print_totuns_fullE()

# ------------------------------------------------------------------
