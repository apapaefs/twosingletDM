from sys import argv
import random
import math
import os
from sys import argv
from generate_trsm_info import * # TRSM info generator (branching ratios, mixing matrices, etc.)
from test_trsm_evolution import * # RGE evolution
#from twosinglet_sigmahhh_read_pickle_nn import * # Triple Higgs Cross Section
from test_trsm_higgstools import * # HiggsTools setup 
from trsm_kstoalphas import * # Convert from k1, k2, k3 to a12, a13, a23
from generate_mg5_trsm_xsecs import * # call MG5 to get the cross section for a specific proces. Make sure that the process has been generated (and check run card for energy/cuts etc!)
from test_trsm_theory_constraints import * # unitarity/boundedness from below
from prettytable import PrettyTable
from datetime import date
from singlet_EWPO import * # Electroweak Precision Observables

###########################################################
# Handle the input here.
# Random seed is the only input
###########################################################

if len(sys.argv) < 1:
    print('generate_trsm_points.py [seed]')
    exit()

ini_seed=int(argv[1])

#############
# OPTIONS
##############

# print debug?
debug = True

# run MG5 on points that pass constraints?
RunMG5 = True

# for random scan within ranges, how many points to run
nrandom=1000

# Energy for xsec calculation
Energy = 13.6

# HH SM XSEC at Energy [pb] (nn23lo1 pdf in MadGraph)
xsec_sm = {}
xsec_sm[13] = 0.01452
xsec_sm[13.6] = 0.01617

# set the Higgs mass:
mhiggs = 125.

# TAG for RUN output
RunTag = str(Energy) + '-' + str(date.today()).replace('-','') + '-' + str(ini_seed) + '-' + str(RunMG5)

# Directory for output:
OutputDir = 'output/'

# reset the output before next run?
ResetOutput = True

# Print additional TRSM point info?
PRINTINFO = True

###########################################################
# some functions
###########################################################

# print the parameter point info info:
def print_info(vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3):
    
    tbl = PrettyTable(["var", "value"])
    tbl.add_row(['vs', vs])
    tbl.add_row(['vx', vx])
    tbl.add_row(['a12', a12])
    tbl.add_row(['a13', a13])
    tbl.add_row(['a23', a23])
    tbl.add_row(['c1', math.cos(a12)])
    tbl.add_row(['c2', math.cos(a13)])
    tbl.add_row(['c3', math.cos(a23)])
    tbl.add_row(['s1', math.sin(a12)])
    tbl.add_row(['s2', math.sin(a13)])
    tbl.add_row(['s3', math.sin(a23)])
    tbl.add_row(['M2', M2])
    tbl.add_row(['w2', w2])
    tbl.add_row(['M3', M3])
    tbl.add_row(['w3', w3])
    tbl.add_row(['k1', k1])
    tbl.add_row(['k2', k2])
    tbl.add_row(['k3', k3])
    tbl.add_row(['K111', K111])
    tbl.add_row(['K112', K112])
    tbl.add_row(['K113', K113])
    tbl.add_row(['K123', K123])
    tbl.add_row(['K122', K122])
    tbl.add_row(['K1111', K1111])
    tbl.add_row(['K1112', K1112])
    tbl.add_row(['K1113', K1113])
    tbl.add_row(['K133', K133])
    print(tbl)
    #print('\n')

def print_constraints(evo, thc, hb, hs, ewpo, wmass):
    tbl = PrettyTable(["Constraint", "Pass/Fail"])
    constraint = {}
    if evo == True:
        constraint['evo'] = 'Pass'
    else:
        constraint['evo'] = 'Fail'
    if thc == True:
        constraint['thc'] = 'Pass'
    else:
        constraint['thc'] = 'Fail'
    if hb == True:
        constraint['hb'] = 'Pass'
    else:
        constraint['hb'] = 'Fail'
    if hs == True:
        constraint['hs'] = 'Pass'
    else:
        constraint['hs'] = 'Fail'
    if ewpo == True:
        constraint['ewpo'] = 'Pass'
    else:
        constraint['ewpo'] = 'Fail'
    if wmass == True:
        constraint['wmass'] = 'Pass'
    else:
        constraint['wmass'] = 'Fail'
    for key in constraint.keys():
        tbl.add_row([key, constraint[key]])
    print(tbl)

# write the point
def write_valid_point_xsec(runtag, m2, m3, vs, vx, a12, a13, a23, xsec, resfrac):
    outfile = OutputDir + 'trsm_points_' + runtag + '.dat'
    filestream = open(outfile,'a')
    filestream.write(str(m2) + '\t' + str(m3) + '\t' + str(vs) + '\t' + str(vx) + '\t' + str(a12) + '\t' + str(a13) + '\t' + str(a23) + '\t' + str(xsec) + '\t' + str(resfrac) + '\n')
    filestream.close()

# write the point without the xsec
def write_valid_point(runtag, m2, m3, vs, vx, a12, a13, a23):
    outfile = OutputDir + 'trsm_points_' + runtag + '.dat'
    filestream = open(outfile,'a')
    filestream.write(str(m2) + '\t' + str(m3) + '\t' + str(vs) + '\t' + str(vx) + '\t' + str(a12) + '\t' + str(a13) + '\t' + str(a23) + '\n')
    filestream.close()

def reset_output(runtag):
    outfile = OutputDir + 'trsm_points_' + runtag + '.dat'
    filestream = open(outfile,'w')
    filestream.close()

    
# MAIN FUNCTION:
def evaluate_trsm_point(myseed, m2_val, m3_val, vs_val, vx_val, a12, a13, a23, runmg5=False):
    # get the point information (widths, scalar couplings)
    vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3, h1_BRs, h2_BRs, h3_BRs, xs136_lo_h1, xs136_lo_h2, xs136_lo_h3 = generate_lams(myseed, m2_val, m3_val, vs_val, vx_val, a12, a13, a23, PRINTINFO)
    Lambdas =[K111,K112,K113,K123,K122,K1111,K1112,K1113,K133]
    if debug is True:
        print_info(vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3)

    # check EWPO
    sinth = np.sin(a12)
    #EWPO_cur = check_EWPO(125.09, M2, sinth, Mz, Mw, Delta_S_central, Delta_T_central, errS, errT, covST) #current # U=0
    EWPO_cur = check_EWPO_wU(125.09, M2, sinth, Mz, Mw, Delta_S_central_wU, Delta_T_central_wU, Delta_U_central_wU, errS_wU, errT_wU, errU_wU, covST_wU, covSU_wU, covTU_wU) #current # U!=0

    #if EWPO_cur != EWPO_cur_wU:
    #    print("EWPOWARNING", EWPO_cur, EWPO_cur_wU)
    

    # check W mass:
    wmass = check_wmass_tania(M2, sinth)
    
        
    # check HiggsTools:
    hb, hs = analyze_parampoint(pred, H1, H2, H3, 125.09, M2, M3, k1, k2, k3, h1_BRs, h2_BRs, h3_BRs)
    if debug is False:
        if hb is False or hs is False:
            return 0
    #thc = theory_constraints(vs, vx, M2, M3, a12, a13, a23)
    #if debug is False:
    #    if thc is False:
    #        return 0
    # test the cosmological constraints
    #evo = test_evo(vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133)
    #if debug is False:
    #    if evo is False:
    #        return 0

    # temporarily removing theory constraints: to be REINSTATED!
    evo = True
    thc = True
    if debug is True:
        print_constraints(evo, thc, hb, hs, EWPO_cur, wmass)
    # get the hh cross section
    # if all constraints are ok, check the xsec for hhh:
    if evo is True and thc is True and hb is True and hs is True and EWPO_cur is True and wmass is True:
        if debug is False:
            print_info(vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3)
            print_constraints(evo, thc, hb, hs)
        if runmg5 is True:
            print('All constraints passed, running mg5 to test cross section, please wait!')
            mg5xsec = get_mg5_xsec('hhh', 'SCAN' + str(Energy), Lambdas, k1, k2, k3, M2, w2, M3, w3,ecm=Energy)
            print('MG5 hh xsec [pb] =', mg5xsec)
            factorxsec = xs136_lo_h2 * h2_BRs[11]
            resfrac = factorxsec/mg5xsec
            print('Factorized h2>h1h1 [pb] =', factorxsec)
            print('Resonant fraction = ', resfrac)
            write_valid_point_xsec(RunTag, m2_val, m3_val, vs_val, vx_val, a12, a13, a23, mg5xsec/xsec_sm[Energy],resfrac)

        else:
            write_valid_point(RunTag, m2_val, m3_val, vs_val, vx_val, a12, a13, a23)
        return 1
    return 0

# round to sgf significant figures
def round_signif(m2, m3, vs, vx, a12, a13, a23, sgf):
    return round_sig(m2, sgf), round_sig(m3, sgf), round_sig(vs, sgf), 0, round_sig(a12, sgf), 0, 0


# random number either -1 or 1:
def randsign():
    return 1 if random.random() < 0.5 else -1


############################################################
# define ranges here
############################################################

# ranges of masses m2 and m3 to scan over
m2_min=255
m2_max=450

m3_min=350
m3_max=600

# ranges of vevs
vs_min=0
vs_max=1000

# no range for vx (DM candidate)!
vx_min=0
vx_max=0

# ranges of k1, k2, k3 (can be positive or negative)
k1_min = 0.95
k1_max = 1.0

# for the grid scan:
num_m2 = 2
num_m3 = 2


############################################################
# Scan begins here
############################################################

print('\nScanning TRSM parameter space')
# reset the output file?
if ResetOutput is True:
    reset_output(RunTag)


# SCAN:
passcounter = 0 # count number of passing points
print('Generating points randomly within ranges')
random.seed(ini_seed)
for i in tqdm(range(0,nrandom)):
        k1=random.uniform(k1_min,k1_max)
        k2=np.sqrt(1-k1**2)
        k3=0
        vx=random.uniform(vx_min, vx_max)
        vs=random.uniform(vs_min, vs_max)
        m2=random.uniform(m2_min, m2_max)
        #if m2 > m3_min:
        #    m3_lim = m2
        #else:
        #    m3_lim = m3_min
        m3=random.uniform(m2+mhiggs, m3_max)
        #a12, a13, a23 = ks_to_angles(k1,k2,k3) 
        a12 = np.arccos(k1)
        a13 = 0
        a23 = 0


        # round to 4 significant figures:
        m2, m3, vs, vx, a12, a13, a23 = round_signif(m2, m3, vs, vx, a12, a13, a23, 4)

        # DM X: set vx = 1E-10 so that we avoid division by 0:
        vx = 1E-10
        
        
        # evaluate: evalpoint is 1 if point passes, 0 if not
        evalpoint = evaluate_trsm_point(ini_seed, m2, m3, vs, vx, a12, a13, a23,runmg5=RunMG5)
        # count the passing points:
        passcounter = passcounter + evalpoint
print('Generated', nrandom,'points, out of which', passcounter, 'are viable')



