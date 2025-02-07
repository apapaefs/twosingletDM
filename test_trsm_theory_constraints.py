import math
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import numpy as np
from math import *
from scipy import interpolate

# check Unitarity/boundedness 
def theory_constraints(vs, vx, M2, M3, a12, a13, a23):

    vh=246
    m1=125.09
    m2=M2
    m3=M3
    
    R = [[0. for x in range(3)] for y in range(3)]
    #calculate the cos functions
    c1 = math.cos(a12)
    c2 = math.cos(a13)
    c3 = math.cos(a23)
    # get the sin functions
    s1 = math.sin(a12)
    s2 = math.sin(a13)
    s3 = math.sin(a23)

    # calculate the elements
    R[0][0] = c1 * c2
    R[0][1] = - s1 * c2
    R[0][2] = - s2
    R[1][0] = s1 * c3 - c1 * s2 * s3
    R[1][1] = c1 * c3 + s1 * s2 * s3
    R[1][2] = - c2 * s3
    R[2][0] = c1 * s2 * c3 + s1 * s3
    R[2][1] = c1 * s3 - s1 * s2 * c3
    R[2][2] = c2 * c3

    myPi=math.pi

    lPhi  =(1/(2*vh**2))*(m1**2*R[0][0]**2+m2**2*R[1][0]**2+m3**2*R[2][0]**2)
    lS =(1/(2*vs**2))*(m1**2*R[0][1]**2+m2**2*R[1][1]**2+m3**2*R[2][1]**2)
    lX =(1/(2*vx**2))*(m1**2*R[0][2]**2+m2**2*R[1][2]**2+m3**2*R[2][2]**2)

    lPhiS =(1/(vh*vs))*(m1**2*R[0][0]*R[0][1]+m2**2*R[1][0]*R[1][1]+m3**2*R[2][0]*R[2][1])
    lPhiX =(1/(vh*vx))*(m1**2*R[0][0]*R[0][2]+m2**2*R[1][0]*R[1][2]+m3**2*R[2][0]*R[2][2])
    lSX =(1/(vs*vx))*(m1**2*R[0][1]*R[0][2]+m2**2*R[1][1]*R[1][2]+m3**2*R[2][1]*R[2][2])

    coeff_x_3=1
    coeff_x_2=-12*lPhi-6*lS-6*lX
    coeff_x=72*lPhi*(lS+lX)-4*(lPhiS**2+lPhiX**2) + 36*lS*lX - lSX**2
    coeff_0=12*lPhi*lSX**2 +24*lPhiS**2*lX + 24*lPhiX**2*lS - 8*lPhiS*lPhiX*lSX -432*lPhi*lS*lX
    coeff_array=[coeff_x_3, coeff_x_2, coeff_x, coeff_0]
    a_array=np.roots(coeff_array)
    success_flag=False
    if(lPhi< 4*myPi and lPhiS< 8*myPi and lPhiX< 8*myPi and lSX<8*myPi):
        if(a_array[0]<16*myPi and a_array[1]<16*myPi and a_array[2]<16*myPi):
            success_flag=True
    return success_flag         







