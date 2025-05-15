#! /usr/bin/env python
import math
import numpy as np
from scipy.interpolate import interp1d

# first some EW parameters:
mz = 91.1876 # from the PDG: http://pdg.lbl.gov/2019/tables/rpp2019-sum-gauge-higgs-bosons.pdf
mw = 80.385
Gf = 1.1663787E-5
alpha = 7.2973525693E-3
v0 = 246.
# the SM value of the self-coupling and Higgs boson mass in [GeV]
v0 = 246.
MH = 125.0
mt = 172.44 # PDG as of 2019
Mw = 80.385
Mz = 91.1876
Mt = 172.44
Gf = 1.1663787E-5;
v0 = math.sqrt(1/math.sqrt(2.)*(1/Gf));
yt = mt/v0 * math.sqrt(2.);
g0sq = 4. * math.sqrt(2.) * Gf * Mw**2;
g1sq = g0sq * (Mz**2 - Mw**2)/Mw**2
g1 = math.sqrt(g1sq)
g2 = math.sqrt(g0sq)
alpha = 7.2973525693E-3

################################################
# W MASS CONSTRAINT                            #
################################################
# the W mass constraint (Tania Robens, Snowmass white paper)
datalist = [ ( np.loadtxt('datafiles/Tania_MW_SnowmassWhitepaper.dat'), 'total' ) ]
for data, label in datalist:
        x1w = data[:,0]
        y1w = data[:,1]
# interpolate the W mass constraint
interpkind = "cubic"
sinthetalim_wmass = interp1d(x1w, y1w, kind=interpkind, bounds_error=False)

def check_wmass_tania(mh2, sinth):
    sinthlim = sinthetalim_wmass(mh2)
    if sinth > sinthlim:
        return False
    else:
        return True


# define the functions for calculation the changes S and T parameters:
# see hep-ph/9409380 appendix C
# x = mh**2 / mz**2
def B(x):
    resB = -1.
    if x >= 4:
        return math.sqrt(x * (x-4)) * math.log( 2/(math.sqrt(x) + math.sqrt(x-4)) )
    elif x < 4 and x >= 0:
        return math.sqrt(x * (4-x)) * math.atan( math.sqrt(4/x -1 ) )
    elif x < 0:
        print("error, B(x) evaluated at x<0")
        exit()

def H_S(x):
    return (3./8.) * x - (1./12.) * x**2 + ( (3.-x)/4. + x**2/24. + 3./(4*(1-x)) ) * x * math.log(x) + (1 - x/3. + x**2/12.) * B(x)

# c**2 = mw**2 / mz**2 
def H_T(x,c):
    return (3. * x / 4.) * ( math.log(x) / (1.-x) - math.log(x/c**2)/(1.-x/c**2) )

# c**2 = mw**2 / mz**2  
def H_U(x,c):
    return - H_S(x) + H_S(x/c**2)


# the "observables"
def O_S(x):
    return (1./math.pi)* H_S(x)

def O_T(x, c, expansion_param):
    return expansion_param * H_T(x, c)

def O_U(x):
    return (1./math.pi)* H_U(x)

def Delta_S(m1, m2, sintheta, mz, mw):
    x1 = m1**2/mz**2
    x2 = m2**2/mz**2
    return sintheta**2 * (O_S(x2) - O_S(x1))

def Delta_T(m1, m2, sintheta, mz, mw):
    expansion_param = Gf * mz**2 / (2 * math.sqrt(2) * math.pi**2 * alpha )
    c = mw/mz
    x1 = m1**2/mz**2
    x2 = m2**2/mz**2
    return sintheta**2 * (O_T(x2, c, expansion_param) - O_T(x1,c, expansion_param))

def Delta_U(m1, m2, sintheta, mz, mw):
    x1 = m1**2/mz**2
    x2 = m2**2/mz**2
    return sintheta**2 * (O_U(x2) - O_U(x1))

# calculate the chisq for two correlated distributions, given the covariance off-diagonal elements cov12
def calc_chisquared_correlated(delta1, delta2, err1, err2, cov12, delta_central1, delta_central2):
    # calculate the inverse of the matrix V required for the chisq:
    Vinv = [[0 for xx in range(2)] for yy in range(2)]
    # now get the chisq:
    detfac = 1./(err1**2 * err2**2)/(1 - cov12**2)
    Vinv[0][0] = detfac * err2**2
    Vinv[0][1] = - detfac * cov12 * err1 * err2
    Vinv[1][0] = - detfac * cov12 * err1 * err2
    Vinv[1][1] = detfac * err1**2
    J0 = (delta_central1-delta1)
    J1 = (delta_central2-delta2)
    chisq_sum = J0 * Vinv[0][0] * J0 + J0 * Vinv[0][1] * J1 + J1 * Vinv[1][1] * J1 + J1 * Vinv[1][0] * J0
    return chisq_sum

# function to get the total chi_squared given m2, m1, sintheta, mz, mw, S, T central, S, T errors, covariance:
def get_chisq_EWPO(m1, m2, sintheta, mz, mw, Sc, Tc, errS, errT, covST):
    deltaS = Delta_S(m1, m2, sintheta, mz, mw)
    deltaT = Delta_T(m1, m2, sintheta, mz, mw)
    #print 'deltaS, deltaT=', deltaS, deltaT
    return calc_chisquared_correlated(deltaS, deltaT, errS, errT, covST, Sc, Tc)

# function to check whether the chi-sq from EWPO excludes or not  (at 2sigma)
def check_EWPO(m1, m2, sintheta, mz, mw, Sc, Tc, errS, errT, covST):
    chisq = get_chisq_EWPO(m1, m2, sintheta, mz, mw, Sc, Tc, errS, errT, covST)
    if chisq > 5.99:
        return False
    else:
        return True

# calculate the chisq for two correlated distributions, given the covariance off-diagonal elements cov12
def calc_chisquared_correlated_wU(delta1, delta2, delta3, err1, err2, err3, cov12, cov13, cov23, delta_central1, delta_central2, delta_central3):
    # see https://arxiv.org/pdf/1407.5342 eq. 17 onwards:
    sigma = np.array([err1, err2, err3]) # sigma_i
    rho = np.array( [[1, cov12, cov13], [cov12, 1, cov23], [cov13, cov23, 1]] ) # rho_ij
    sigmasq = np.dot(sigma,np.dot(rho,sigma)) # sigma_ij^2
    sigmasqInv = np.linalg.inverse(sigmasq) # (sigma_ij^2)^-1
    deltaOmOc = np.array([(delta_central1-delta1), (delta_central2-delta2), (delta_central3-delta3)]) # DeltaO - DeltaO_central
    # delta_chisq: 
    chisq_sum = np.dot( deltaOmOc, np.dot(sigmasqInv, deltaOmOc))
    return chisq_sum

# function to get the total chi_squared given m2, m1, sintheta, mz, mw, S, T central, S, T errors, covariance:
# with U!=0
def get_chisq_EWPO_wU(m1, m2, sintheta, mz, mw, Sc, Tc, Uc, errS, errT, errU, covST, covSU, covTU):
    deltaS = Delta_S(m1, m2, sintheta, mz, mw)
    deltaT = Delta_T(m1, m2, sintheta, mz, mw)
    deltaU = Delta_U(m1, m2, sintheta, mz, mw)
    return calc_chisquared_correlated_wU(deltaS, deltaT, deltaU, errS, errT, errU, covST, covSU, covTU, Sc, Tc, Uc)


# function to check whether the chi-sq from EWPO excludes or not  (at 2sigma)
def check_EWPO_wU(m1, m2, sintheta, mz, mw, Sc, Tc, Uc, errS, errT, errU, covST, covSU, covTU):
    chisq = get_chisq_EWPO_wU(m1, m2, sintheta, mz, mw, Sc, Tc, errS, errT, covST)
    if chisq > 7.82: # three degrees of freedom! 
        return False
    else:
        return True


   
# the GFITTER CURRENT central values and errors, see https://arxiv.org/pdf/1407.3792.pdf, page 11

# for U = 0:

# central values:
Delta_S_central = 0.06
Delta_T_central = 0.10
# errors:
errS = 0.09
errT = 0.07
# correlation (rho_12):
covST=0.91
# the GFITTER FUTURE central values and errors, see https://arxiv.org/pdf/1407.3792.pdf, page 13, Table 3
# central values:
Delta_S_central_F = 0.00
Delta_T_central_F = 0.00
# errors:
errS_F = math.sqrt(0.017**2 + 0.006**2) # experimental + theoretical errors in quadrature
errT_F = math.sqrt(0.022**2 + 0.005**2) 
# correlation (rho_12):
covST_F=0.91


# for U!=0:

# central values:
Delta_S_central_wU = 0.05
Delta_T_central_wU = 0.09
Delta_U_central_wU = 0.01
# errors:
errS_wU = 0.11
errT_wU = 0.13
errU_wU = 0.11
# correlation:
covST_wU=0.90
covSU_wU=-0.59
covTU_wU=-0.83

    
test = False
if test is True:
    import mpmath as mp
    import numpy as np
    import pylab as pl
    from scipy import interpolate, signal
    import matplotlib.font_manager as fm
    from datetime import date
    from matplotlib.ticker import MultipleLocator
    import matplotlib.patches as mpatches
    import scipy.interpolate
    import matplotlib.gridspec as gridspec
    import matplotlib
    import pylab as pl
    import matplotlib.pyplot as plt

    #########################
    # test the whole thing: #
    #########################
    # EW parameters:
    # G_F * m_z**2 / (2 * math.sqrt(2) * math.pi**2 * alpha) = 0.4761 from hep-ph/9409380
    #expansion_parameter = 0.4761 [compare to hep-ph/9409380]
    
    
    # the GFITTER CURRENT central values and errors, see https://arxiv.org/pdf/1407.3792.pdf, page 11
    
    # for U = 0:
    # central values:
    Delta_S_central = 0.06
    Delta_T_central = 0.10
    # errors:
    errS = 0.09
    errT = 0.07
    # correlation (rho_12):
    covST=0.91

    
    
    # the GFITTER FUTURE central values and errors, see https://arxiv.org/pdf/1407.3792.pdf, page 13, Table 3
    # central values:
    Delta_S_central_F = 0.00
    Delta_T_central_F = 0.00
    # errors:
    errS_F = math.sqrt(0.017**2 + 0.006**2) # experimental + theoretical errors in quadrature
    errT_F = math.sqrt(0.022**2 + 0.005**2) 
    # correlation (rho_12):
    covST_F=0.91



    # get the values predicted for specific m2, m1, sintheta

    print('expansion param.=', Gf * mz**2 / (2 * math.sqrt(2) * math.pi**2 * alpha ))

    # now pick a Benchmark point:
    #  | name |  #         |   mh1   |  mh2   |    Gh1    |   Gh2    |  stheta  |  ctheta | BR(h2->h1h1) | HB res. |  HS res.  | XS13(mh2)[pb] | XS14(mh2)[pb] | XS8(mh2)[pb] | XS7(mh2)[pb] |
    #   B4min |  9         | 124.91 | 463.13 | 0.0040673 | 0.093934 | 0.02576  | 0.99967 |    0.5825    |    1    |  0.967328 |   0.0013471   |   0.0045521   |  0.0011347   |  0.00078302  |

    m1 = 124.91
    m2 = 463.13
    sintheta = 0.02576
    costheta = math.sqrt(1-sintheta**2)
    print('chi-sq for', m1, m2, sintheta, costheta, '=', (m1, m2, sintheta, mz, mw, Delta_S_central, Delta_T_central, errS, errT, covST))
    print('chi-sq (future) for', m1, m2, sintheta, costheta, '=', get_chisq_EWPO(m1, m2, sintheta, mz, mw, Delta_S_central_F, Delta_T_central_F, errS_F, errT_F, covST_F))


    # scan and create exclusion plot:
    # Table 39.2 of http://pdg.lbl.gov/2019/reviews/rpp2019-rev-statistics.pdf:
    # require Delta Chi^2 = 2.30, 6.18 (5.99 for 95% exactly), 11.83 for 1sigma, 2sigma, 3sigma confidence-level.

    # fix m1 to = 125 GeV
    m1 = 125.
    m2_array = np.arange(150, 1000, 5) 
    costheta_array = np.arange(0.8, 1.0, 0.001)
    # arrays for the plot
    m2ar = []
    ctar = []
    chisqar = []
    chisqar_F = []
    # loop over and calculate chisq
    for m2i in range(len(m2_array)):
        for cthi in range(len(costheta_array)):
            sthi = math.sqrt(1-costheta_array[cthi]**2)
            chisqi = get_chisq_EWPO(m1, m2_array[m2i], sthi, mz, mw, Delta_S_central, Delta_T_central, errS, errT, covST) # present
            chisqi_F  = get_chisq_EWPO(m1, m2_array[m2i], sthi, mz, mw, Delta_S_central_F, Delta_T_central_F, errS_F, errT_F, covST_F) # future
            chisqar.append(chisqi)
            chisqar_F.append(chisqi_F)
            m2ar.append(m2_array[m2i])
            ctar.append(costheta_array[cthi])
            #print m2_array[m2i], costheta_array[cthi], sthi, chisqi
            
    ############        
    # plot     #
    ############
    gs = gridspec.GridSpec(6,6)
    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.grid(False)

    xi = np.arange(150, 1000, 1)
    yi = np.arange(0.8, 1.0, 0.0005)
    zi = matplotlib.mlab.griddata(m2ar, ctar, chisqar, xi, yi, interp='linear')
    #cs = plt.contourf(xi, yi, zi, levels=[0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4.5, 5, 5.5, 6, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11, 11.5, 12])
    cs = plt.contour(xi, yi, zi, levels=[2.30, 6.18 , 11.83], extend='both', colors='k')
    #cbar = fig.colorbar(cs)
    #cbar.ax.get_yaxis().labelpad = 15
    #cbar.ax.set_ylabel('$\\chi^2$', rotation=270)
    #plt.clabel(cs, inline=1, fontsize=10)
    manual_locations = [(300, 0.96), (500, 0.89)]

    strs = ['$2\\sigma$', '$3\\sigma$']
    fmt = {}
    for l, s in zip(cs.levels, strs):
        fmt[l] = s
    ax.clabel(cs, cs.levels, inline=True, fmt=fmt, fontsize=12, manual=manual_locations)
    ax.set_ylabel('$\\cos \\theta$', fontsize=20)
    ax.set_xlabel('$m_2$ [GeV]', fontsize=20)
    ax.set_title('EWPO exclusion (current)')
    ax.set_ylim(0.8,0.99)
    ax.set_xlim(200,950)

    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(50))

    ax.yaxis.set_major_locator(MultipleLocator(0.02))
    ax.yaxis.set_minor_locator(MultipleLocator(0.005))
    
    plot_type = 'EWPO_current'
    outputdirectory = 'plots_210420/'
    # save the figure
    print('saving the figure')
    # save the figure in PDF format
    infile = plot_type + '.dat'
    print('---')
    print('output in', outputdirectory + infile.replace('.dat','.pdf'))
    pl.savefig(outputdirectory + infile.replace('.dat','.pdf'), bbox_inches='tight')
    pl.savefig(outputdirectory + infile.replace('.dat','.png'), bbox_inches='tight', scale=0.1)
    pl.close(fig)

    ##################      
    # plot  FUTURE   #
    ##################
    gs = gridspec.GridSpec(6,6)
    fig = pl.figure()
    ax = fig.add_subplot(111)
    ax.grid(False)

    xi = np.arange(150, 1000, 1)
    yi = np.arange(0.8, 1.0, 0.0005)
    zi = matplotlib.mlab.griddata(m2ar, ctar, chisqar_F, xi, yi, interp='linear')
    #cs = plt.contourf(xi, yi, zi, levels=[0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4.5, 5, 5.5, 6, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0, 10.5, 11, 11.5, 12])
    cs = plt.contour(xi, yi, zi, levels=[2.30, 6.18 , 11.83], extend='both', colors='k')
    #cbar = fig.colorbar(cs)
    #cbar.ax.get_yaxis().labelpad = 15
    #cbar.ax.set_ylabel('$\\chi^2$', rotation=270)
    manual_locations = [(300, 0.96), (350, 0.94), (440, 0.89)]
    strs = ['$1\\sigma$', '$2\\sigma$', '$3\\sigma$']
    fmt = {}
    for l, s in zip(cs.levels, strs):
        fmt[l] = s
    plt.clabel(cs, inline=1, fontsize=12, fmt=fmt,manual=manual_locations)

    ax.set_ylabel('$\\cos \\theta$', fontsize=20)
    ax.set_xlabel('$m_2$ [GeV]', fontsize=20)
    ax.set_title('EWPO exclusion (ILC/GigaZ)')
    ax.set_ylim(0.8,0.99)
    ax.set_xlim(200,950)

    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(50))

    ax.yaxis.set_major_locator(MultipleLocator(0.02))
    ax.yaxis.set_minor_locator(MultipleLocator(0.005))
    
    plot_type = 'EWPO_future'
    outputdirectory = 'plots_210420/'
    # save the figure
    print('saving the figure')
    # save the figure in PDF format
    infile = plot_type + '.dat'
    print('---')
    print('output in', outputdirectory + infile.replace('.dat','.pdf'))
    pl.savefig(outputdirectory + infile.replace('.dat','.pdf'), bbox_inches='tight')
    pl.savefig(outputdirectory + infile.replace('.dat','.png'), bbox_inches='tight', scale=0.1)
    pl.close(fig)

