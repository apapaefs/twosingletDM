import math
from math import log10, floor
import Higgs.predictions as HP
import Higgs.bounds as HB
import Higgs.signals as HS
import numpy as np
import pandas as pd
from prettytable import PrettyTable
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator
from tqdm import tqdm
from matplotlib import ticker, cm
from scipy import interpolate
from scipy.optimize import fsolve, root
from functools import partial

def round_sig(x, sig=2):
    if x == 0.:
        return 0.
    if math.isnan(x) is True:
        print('Warning, NaN!')
        return 0.
    return round(x, sig-int(floor(log10(abs(x))))-1)


##################################
# SET UP HIGGS TOOLS:
#####################################

mhSM = 125.09

pred = HP.Predictions() # create the model predictions

bounds = HB.Bounds('./hbdataset') # load HB dataset
signals = HS.Signals('./hsdataset') # load HS dataset

# add a SM-like Higgs boson with SM-like couplings

H1 = pred.addParticle(HP.BsmParticle("H1", "neutral", "even"))

H1.setMass(mhSM)
HP.effectiveCouplingInput(H1, HP.scaledSMlikeEffCouplings(1.0))
# get the SM chi-squared for HiggsSignals
ress_SM = signals(pred)

print("HiggsSignals chi-sq. for SM=", ress_SM)

# add second BSM Higgs boson 

H2 = pred.addParticle(HP.BsmParticle("H2", "neutral", "even"))

# add third BSM Higgs boson

H3 = pred.addParticle(HP.BsmParticle("H3", "neutral", "even"))


##########################################
# FUNCTIONS 
##########################################

# ensure that the BR arrays sum up to unity to avoid issues with HiggsTools
def ensure_sum_unit(brarray):
    sumarray = brarray.sum() - brarray[-1]
    brarray_new = np.array([None] * len(brarray))
    for i in range(len(brarray)-1):
        brarray_new[i] = round_sig(brarray[i]/sumarray,100)
    brarray_new[-1] = round_sig(brarray[-1],3)
    #print(brarray_new.sum() - brarray_new[-1])
    return brarray_new

def analyze_parampoint(pred, H1, H2, H3, M1, M2, M3, k1, k2, k3, h1_BRs, h2_BRs, h3_BRs):

    h1_BRs = ensure_sum_unit(h1_BRs)
    h2_BRs = ensure_sum_unit(h2_BRs)

    # set the present masses for the scalars
    H1.setMass(M1)
    H2.setMass(M2)
    # get the branching ratios
    HP.effectiveCouplingInput(H1, HP.scaledSMlikeEffCouplings(k1))
    HP.effectiveCouplingInput(H2, HP.scaledSMlikeEffCouplings(k2))
    # RESET BRs BEFORE SETTING THEM TO AVOID ISSUES WITH BR>1
    # H -> bbbar | H -> tautau | H -> mumu | H -> cc | H -> ss | H -> tt | H -> gg | H -> gammagamma | H -> Zgamma | H -> WW | H -> ZZ
    # H1
    H1.setBr('bb', 0.)
    H1.setBr('tautau', 0.)
    H1.setBr('mumu', 0.)
    H1.setBr('cc', 0.)
    H1.setBr('ss', 0.)
    H1.setBr('tt', 0.)
    H1.setBr('gg', 0.)
    H1.setBr('gamgam', 0.)
    H1.setBr('Zgam', 0.)
    H1.setBr('WW', 0.)
    H1.setBr('ZZ', 0.)
    # SET THE BRS
    H1.setBr('bb', h1_BRs[0])
    H1.setBr('tautau', h1_BRs[1])
    H1.setBr('mumu', h1_BRs[2])
    H1.setBr('cc', h1_BRs[3])
    H1.setBr('ss', h1_BRs[4])
    H1.setBr('tt', h1_BRs[5])
    H1.setBr('gg', h1_BRs[6])
    H1.setBr('gamgam', h1_BRs[7])
    H1.setBr('Zgam',h1_BRs[8])
    H1.setBr('WW', h1_BRs[9])
    H1.setBr('ZZ', h1_BRs[10])
    H1.setTotalWidth(h1_BRs[-1])
    # H2
    # RESET BRs BEFORE SETTING THEM TO AVOID ISSUES WITH BR>1
    H2.setBr('bb', 0.)
    H2.setBr('tautau', 0.)
    H2.setBr('mumu', 0.)
    H2.setBr('cc', 0.)
    H2.setBr('ss', 0.)
    H2.setBr('tt', 0.)
    H2.setBr('gg', 0.)
    H2.setBr('gamgam', 0.)
    H2.setBr('Zgam', 0.)
    H2.setBr('WW', 0.)
    H2.setBr('ZZ', 0.)
    # SET THE BRS
    H2.setBr('bb', h2_BRs[0])
    H2.setBr('tautau', h2_BRs[1])
    H2.setBr('mumu', h2_BRs[2])
    H2.setBr('cc', h2_BRs[3])
    H2.setBr('ss', h2_BRs[4])
    H2.setBr('tt', h2_BRs[5])
    H2.setBr('gg', h2_BRs[6])
    H2.setBr('gamgam', h2_BRs[7])
    H2.setBr('Zgam',h2_BRs[8])
    H2.setBr('WW',h2_BRs[9])
    H2.setBr('ZZ',h2_BRs[10])
    H2.setBr('H1','H1',h2_BRs[11])
    H2.setTotalWidth(h2_BRs[-1])

    # get the HiggsBounds result
    resb = bounds(pred)

    # get the HiggsSignals result
    ress = signals(pred)
    #print(signals(pred).appliedLimits)
    #print(ress)
    HS_allowed = False
    if ress - ress_SM < 4.00:
        HS_allowed = True
    else:
        HS_allowed = False

    return resb.allowed, HS_allowed
