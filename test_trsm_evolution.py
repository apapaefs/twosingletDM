import math
from scipy.integrate import solve_ivp
import numpy as np
from math import *

def rhs(mu, v): 
    return [-7/(16.0*mu*pi*pi)*v[0]*v[0]*v[0],#gg3 
            -19/(96.0*mu*pi*pi)*v[1]*v[1]*v[1],#gg2
            41/(96*mu*pi*pi)*v[2]*v[2]*v[2],#gg1
            v[3]*(-17/12.0*v[2]*v[2]-9/4.0*v[1]*v[1]-8*v[0]*v[0]+9/2.0*v[3]*v[3])/(16*pi*pi*mu),#ytt
            (3*v[2]**4 + 9*v[1]**4 + 6*v[2]**2*(v[1]**2 -4*v[4])-72*v[1]**2*v[4]+96*v[4]*(2*v[4] + v[3]**2) + 4*(v[7]**2 + v[8]**2 -12*v[3]**4))/(128*mu*pi**2),#lH
            (4*v[7]**2 + v[9]**2 +36*v[5]**2)/(32*pi**2*mu),#lS1
            (4*v[8]**2 + v[9]**2 +36*v[6]**2)/(32*pi**2*mu), #lS2
            (-3*v[2]**2*v[7]-9*v[1]**2*v[7]+2*v[8]*v[9]+4*v[7]*(3*v[3]**2+6*v[4]+2*v[7]+3*v[5]))/(32*pi**2*mu),#lM1
            (-3*v[2]**2*v[8]-9*v[1]**2*v[8]+2*v[7]*v[9]+4*v[8]*(3*v[3]**2+6*v[4]+2*v[8]+3*v[6]))/(32*pi**2*mu), #lM2
           (2*v[7]*v[8]+v[9]*(2*v[9]+3*(v[5]+v[6])))/(8*pi**2*mu),#lM3
           (2*(v[7]*v[11]+v[8]*v[12])-3*v[10]*(v[2]**2 + 3*v[1]**2-4*(2*v[4]+v[3]**2)))/(32*pi**2*mu),#mm1
           (4*v[7]*v[10]+6*v[5]*v[11] + v[9]*v[12])/(16*pi**2*mu),#mmS1
           (4*v[8]*v[10]+v[9]*v[11]+6*v[6]*v[12])/(16*pi**2*mu)#mmS2
            ]

def solv_eq(low_esc,# low energy scale for evolution 
            high_esc,#high energy scale for evolution 
            part, #retained for backward-compatible call signatures
            ini_conditions,#initial conditions for equations
            test_esc#test scale to check divergence
            ):
    _ = part
    scales = (low_esc, high_esc, test_esc)
    try:
        low_esc, high_esc, test_esc = (float(scale) for scale in scales)
    except (TypeError, ValueError, OverflowError):
        return 0
    if not all(math.isfinite(scale) for scale in (low_esc, high_esc, test_esc)):
        return 0
    if low_esc == high_esc or not (
        min(low_esc, high_esc) <= test_esc <= max(low_esc, high_esc)
    ):
        return 0

    try:
        res = solve_ivp(
            rhs,
            (low_esc, high_esc),
            ini_conditions,
            t_eval=[test_esc],
        )
    except (FloatingPointError, OverflowError, RuntimeError, ValueError):
        return 0

    if not res.success:
        return 0
    values = np.asarray(res.y)
    if (
        values.ndim != 2
        or values.shape[0] != len(ini_conditions)
        or values.shape[1] != 1
    ):
        return 0
    return int(np.all(np.isfinite(values[:, 0])))



def test_evo(vs, vx, M2, M3, a12, a13, a23, width1, width2, width3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133):

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
    #print('R matrix=',R)

    # "rescaling" coefficients
    k1 = R[0][0]
    k2 = R[1][0]
    k3 = R[2][0]

    # put the lamdbas in a single array
    Lambdas =[K111, K112, K113, K123, K122, K1111, K1112, K1113, K133]

    # masses and widths:
    vh=246
    m1=125.09
    m2=M2
    m3=M3
    w2=width2
    w3=width3


    ###Initial conditions Jorinde's program
    ###Fixed initial conditions
    g1init=sqrt(0.1273)
    g2init=sqrt(0.424)
    g3init=1.221
    yttinit=0.96738
    ##Dynamic initial conditions
    lhinit  =(1/(2*vh**2))*(m1**2*R[0][0]**2+m2**2*R[1][0]**2+m3**2*R[2][0]**2)
    lS1init =(1/(2*vs**2))*(m1**2*R[0][1]**2+m2**2*R[1][1]**2+m3**2*R[2][1]**2)
    lS2init =(1/(2*vx**2))*(m1**2*R[0][2]**2+m2**2*R[1][2]**2+m3**2*R[2][2]**2)
    lM1init =(1/(vh*vs))*(m1**2*R[0][0]*R[0][1]+m2**2*R[1][0]*R[1][1]+m3**2*R[2][0]*R[2][1])
    lM2init =(1/(vh*vx))*(m1**2*R[0][0]*R[0][2]+m2**2*R[1][0]*R[1][2]+m3**2*R[2][0]*R[2][2])
    lM3init =(1/(vs*vx))*(m1**2*R[0][1]*R[0][2]+m2**2*R[1][1]*R[1][2]+m3**2*R[2][1]*R[2][2])
    mm1init =-((vh**2)*lhinit  + (vs**2)/2.0*lM1init + (vx**2)/2.0*lM2init)
    mmS1init=-((vs**2)*lS1init + (vh**2/2.0)*lM1init + (vx**2)/2.0*lM3init)
    mmS2init=-((vx**2)*lS2init + (vh**2/2.0)*lM2init + (vs**2)/2.0*lM3init)

    low_esc=91
    high_esc=1000
    part=0.0001
    test_esc=900

    ini_conditions=[                 g3init,               #v[0]#gg3
                                 g2init,  #v[1]#gg2
                                 g1init, #v[2]#gg1
                                 yttinit,             #v[3]#ytt
                                 lhinit, #v[4]#lH
                                 lS1init, #v[5]#lS1
                                 lS2init,  #v[6]#lS2
                                 lM1init,  #v[7]#lM1
                                 lM2init,  #v[8]#lM2
                                 lM3init,   #v[9]#lM3
                                 mm1init, #v[10]#mm1
                                 mmS1init, #v[11]##mmS1
                                 mmS2init]#v[12]##mmS2




    flag=solv_eq(low_esc,# low energy scale for evolution
            high_esc,#high energy scale for evolution
            part, #size of partition for interpolation
            ini_conditions,#initial conditions for equations
            test_esc#test scale to check divergence
            )
    # True is ok (1), False is not (0)
    if flag == 1:
        flag = True
    else:
        flag = False
    return flag

def test_evo_vxzero(vs, M2, M3, a12, lX, lPhiX, lSX, width1, width2, width3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133):

    # for vx = 0:
    a13 = 0
    a23 = 0
    vx = 0
    
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
    #print('R matrix=',R)

    # "rescaling" coefficients
    k1 = R[0][0]
    k2 = R[1][0]
    k3 = R[2][0]

    # put the lamdbas in a single array
    Lambdas =[K111, K112, K113, K123, K122, K1111, K1112, K1113, K133]

    # masses and widths:
    vh=246
    m1=125.09
    m2=M2
    m3=M3
    w2=width2
    w3=width3


    ###Initial conditions Jorinde's program
    ###Fixed initial conditions
    g1init=sqrt(0.1273)
    g2init=sqrt(0.424)
    g3init=1.221
    yttinit=0.96738
    ##Dynamic initial conditions
    lhinit  =(1/(2*vh**2))*(m1**2*R[0][0]**2+m2**2*R[1][0]**2+m3**2*R[2][0]**2)
    lS1init =(1/(2*vs**2))*(m1**2*R[0][1]**2+m2**2*R[1][1]**2+m3**2*R[2][1]**2)
    #lS2init =(1/(2*vx**2))*(m1**2*R[0][2]**2+m2**2*R[1][2]**2+m3**2*R[2][2]**2)
    lS2init = lX
    lM1init =(1/(vh*vs))*(m1**2*R[0][0]*R[0][1]+m2**2*R[1][0]*R[1][1]+m3**2*R[2][0]*R[2][1])
    #lM2init =(1/(vh*vx))*(m1**2*R[0][0]*R[0][2]+m2**2*R[1][0]*R[1][2]+m3**2*R[2][0]*R[2][2])
    lM2init = lPhiX
    #lM3init =(1/(vs*vx))*(m1**2*R[0][1]*R[0][2]+m2**2*R[1][1]*R[1][2]+m3**2*R[2][1]*R[2][2])
    lM3init = lSX

    mm1init =-((vh**2)*lhinit  + (vs**2)/2.0*lM1init + (vx**2)/2.0*lM2init)
    mmS1init=-((vs**2)*lS1init + (vh**2/2.0)*lM1init + (vx**2)/2.0*lM3init)
    mmS2init=-((vx**2)*lS2init + (vh**2/2.0)*lM2init + (vs**2)/2.0*lM3init)

    low_esc=91
    high_esc=1000
    part=0.0001
    test_esc=900

    ini_conditions=[                 g3init,               #v[0]#gg3
                                 g2init,  #v[1]#gg2
                                 g1init, #v[2]#gg1
                                 yttinit,             #v[3]#ytt
                                 lhinit, #v[4]#lH
                                 lS1init, #v[5]#lS1
                                 lS2init,  #v[6]#lS2
                                 lM1init,  #v[7]#lM1
                                 lM2init,  #v[8]#lM2
                                 lM3init,   #v[9]#lM3
                                 mm1init, #v[10]#mm1
                                 mmS1init, #v[11]##mmS1
                                 mmS2init]#v[12]##mmS2




    flag=solv_eq(low_esc,# low energy scale for evolution
            high_esc,#high energy scale for evolution
            part, #size of partition for interpolation
            ini_conditions,#initial conditions for equations
            test_esc#test scale to check divergence
            )
    # True is ok (1), False is not (0)
    if flag == 1:
        flag = True
    else:
        flag = False
    return flag
