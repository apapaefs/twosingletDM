from scipy.optimize import fsolve
import math

def eq_ks_to_angles(p, k1, k2, k3):
    a12, a13, a23 = p
    c1 = math.cos(a12)
    c2 = math.cos(a13)
    c3 = math.cos(a23)
    s1 = math.sin(a12)
    s2 = math.sin(a13)
    s3 = math.sin(a23)
    return (c1 * c2 - k1, s1 * c3 - c1 * s2 * s3 - k2, c1 * s2 * c3 + s1 * s3 - k3)

def angles_to_ks(a12, a13, a23):
    c1 = math.cos(a12)
    c2 = math.cos(a13)
    c3 = math.cos(a23)
    s1 = math.sin(a12)
    s2 = math.sin(a13)
    s3 = math.sin(a23)
    k1 = c1 * c2
    k2 = s1 * c3 - c1 * s2 * s3
    k3 = c1 * s2 * c3 + s1 * s3
    return k1, k2, k3

def convert_twopi(a12, a13, a23):
    a12 = a12%(2*math.pi)
    a13 = a13%(2*math.pi)
    a23 = a23%(2*math.pi)
    return a12, a13, a23

def ks_to_angles(k1,k2,k3):
    a12, a13, a23 = fsolve(eq_ks_to_angles, (0.5, 0.5, 0.5), args=(k1,k2,k3))
    a12, a13, a23 = convert_twopi(a12, a13, a23)
    return a12, a13, a23


# TEST: 

#k1 = 0.9968128327470546
#k2 = -0.07195383520701071
#k3 = -0.0344502840306862  

#a12, a13, a23 = ks_to_angles(k1,k2,k3)

#print('sol a12, a13, a23=',a12, a13, a23)
#c1 = math.cos(a12)
#c2 = math.cos(a13)
#c3 = math.cos(a23)
#s1 = math.sin(a12)
#s2 = math.sin(a13)
#s3 = math.sin(a23)
#print('sol c1, c2, c3=',c1, c2, c3)
#print('sol s1, s2, s3=',s1, s2, s3)
#print(eq_ks_to_angles((a12, a13, a23), k1, k2, k3))
#k1, k2, k3 = angles_to_ks(a12, a13, a23)
#print('sol k1, k2, k3=', k1, k2, k3)
#a12=0.041126274366412405
#a13=0.06847625639900101
#a23=2.5587754659637483
#k1=0.9968128327470546
#k2=-0.07195383520701071
#k3=-0.0344502840306862
#print('tru a12, a13, a23=', a12, a13, a23)
#c1 = math.cos(a12)
#c2 = math.cos(a13)
#c3 = math.cos(a23)
#s1 = math.sin(a12)
#s2 = math.sin(a13)
#s3 = math.sin(a23)
#print('tru c1, c2, c3=',c1, c2, c3)
#print('tru s1, s2, s3=',s1, s2, s3)
#print('true k1, k2, k3=', k1, k2, k3)

