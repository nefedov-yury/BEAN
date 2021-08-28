#!/usr/bin/python3
# average of two periods of data taking

import numpy as np
#  from scipy.linalg import inv

def test():
    A = np.array([[3, 1], [1, 4]])
    print('A=\n',A)
    # Invers matrix
    invA = np.linalg.inv(A)
    print('invA=\n',invA)
    # check
    print('A*invA=\n',np.dot(A,invA))

# phi-eta: 06jul2021
#   Br(2009)_n = (9.151 ± 0.320 ± 0.242) * 10^{-4}
#   Br(2012)_n = (9.105 ± 0.204 ± 0.187) * 10^{-4}
#
#   Br(2009)_p = (7.944 ± 0.312 ± 0.157) * 10^{-4}
#   Br(2012)_p = (7.904 ± 0.226 ± 0.200) * 10^{-4}
#
#   δBr_com_n = 1.277%
#   δBr_com_p = 1.315%

# KK-eta: 26-30июн21
#   Br(KKη 2009) = (4.556 ± 0.150 ± 0.067) * 10^{-4}
#   Br(KKη 2012) = (4.533 ± 0.086 ± 0.047) * 10^{-4}
#
#   δBr_com = 0.747%

#  debug=True
debug=False

#  Br = np.array([ 9.151, 9.105])
#  Br = np.array([ 7.944, 7.904])
#
Br = np.array([ 4.556, 4.533])
print(' Br= ',Br,'\n')

# correlated
#  sBr_common = 0.01277
#  sBr_common = 0.01315
#
sBr_common = 0.00747
s1 = Br[0]*sBr_common
s2 = Br[1]*sBr_common
Ec = np.array([[s1**2, s1*s2], [s1*s2, s2**2]])
if debug: print(' correlated: Ec=\n',Ec)

# uncorrelated part
#  Eu = np.array([[0.320**2+0.242**2-s1**2, 0], [0, 0.204**2+0.187**2-s2**2]])
#  Eu = np.array([[0.312**2+0.157**2-s1**2, 0], [0, 0.226**2+0.200**2-s2**2]])
#
Eu = np.array([[0.150**2+0.067**2-s1**2, 0], [0, 0.086**2+0.047**2-s2**2]])

if debug:  print(' uncorrelated: Eu=\n',Eu)

E = Eu+Ec
print(' Error Matrix: E=\n',E,'\n')

invE = np.linalg.inv(E)
if debug: print(' invE=\n',invE,'\n')

print(' U is a vector whose components are all unity\n')
U = np.ones((2))
if debug: print('U=\n',U,'\n U.dtype=',U.dtype)

Alpha = np.dot(U,invE)
Alpha = Alpha / np.dot(U.T,Alpha)
print(' Alpha = Einv U /[U(tr) Einv U]\n',Alpha,'\n')

avBr = np.dot(Alpha,Br.T)
print(' <Br> = Alpha(tr) * Br =', avBr, '\n')

sig2 = np.dot(Alpha.T,np.dot(E,Alpha))
print(' <sigma> = sqrt(Alpha(tr) * E * Alpha) =', np.sqrt(sig2) ,'\n')


