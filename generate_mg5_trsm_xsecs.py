import subprocess
import os.path
import math
import numpy as np
from math import log10, floor

# MG5/aMC sub-dir (INCLUDE THE SLASH AT THE END!):
MGLocation = '/home/apapaefs/Projects/TwoSingletDM/twosingletDM/MG5_aMC_v2_9_23/'


# Process sub-dirs (INCLUDE THE SLASH AT THE END!):
ProcLocation = {}
ProcLocation['hh'] = 'gg_hh_twoscalar/'
ProcLocation['hhh'] = 'gg_hhh_twoscalar/'


def round_sig(x, sig=2):
    if x == 0.:
        return 0.
    if math.isnan(x) is True:
        print('Warning, NaN!', x)
        return 0.
    return round(x, sig-int(floor(log10(abs(x))))-1)

# function to run MG5:
def drive_mg(process, runnum, mgloc, k1choice, k2choice, k3choice, LambdasArray, m2, w2, m3, w3, nevents, nruns,output=False,ecm=13):
    if process in ProcLocation:
        procloc = ProcLocation[process]
    else:
        print('Process', process,'is not defined, exiting!')
    filename = mgloc + procloc + '/gg_' + process + '_lambdavar_run' + str(runnum) + '.dcmd'
    #print('generating mg5input:', filename)
    ebeam1 = ecm*1000/2
    ebeam2 = ebeam1
    counter = 0
    for lams in LambdasArray:
        if counter > nruns:
            break
        #print(lams)
        lhe = 'run' + str(runnum) + '_m2_' + str(m2) + '_m3_' + str(m3) + '_w2_' + str(w2) + '_w3_' + str(w3) + '_' + '_'.join((lams)) + '/unweighted_events.lhe.gz'
        lhefile = MGLocation + procloc + 'Events/' + lhe
        #print('lhefile=', lhefile)
        #TestBool = True
        #if TestBool is False:
        while os.path.exists(lhefile) is False:
            filestream = open(filename,'w')
            #filestream.write('launch run' + str(RunNum) + '_m2_' + str(m2) + '_m3_' + str(m3) + '_w2_' + str(w2) + '_w3_' + str(w3) + '_' + '_'.join((lams)) + ' --accuracy=0.25 --points=300 --iterations=1\n0\n')
            filestream.write('generate_events run' + str(runnum) + '_m2_' + str(m2) + '_m3_' + str(m3) + '_w2_' + str(w2) + '_w3_' + str(w3) + '_' + '_'.join((lams)) + ' --accuracy=0.25 --points=300 --iterations=1\n')
            filestream.write('set ebeam1 ' + str(ebeam1) + '\n')
            filestream.write('set ebeam2 ' + str(ebeam2) + '\n')
            filestream.write('set Meta ' + str(m2) + '\n')
            filestream.write('set Weta ' + str(w2) + '\n')
            filestream.write('set Miota ' + str(m3) + '\n')
            filestream.write('set Wiota ' + str(w3) + '\n')
            filestream.write('set k1 ' + str(k1choice) + '\n')
            filestream.write('set k2 ' + str(k2choice) + '\n')
            filestream.write('set k3 ' + str(k3choice) + '\n')
            filestream.write('set kap111 ' + str(lams[0]) + '\n')
            filestream.write('set kap112 ' + str(lams[1]) + '\n')
            filestream.write('set kap113 ' + str(lams[2]) + '\n')
            filestream.write('set kap123 ' + str(lams[3]) + '\n')
            filestream.write('set kap122 ' + str(lams[4]) + '\n')
            filestream.write('set kap1111 ' + str(lams[5]) + '\n')
            filestream.write('set kap1112 ' + str(lams[6]) + '\n')
            filestream.write('set kap1113 ' + str(lams[7]) + '\n')
            filestream.write('set kap133 ' + str(lams[8]) + '\n')
            filestream.write('set nevents ' + str(nevents) + '\n')
            filestream.write('0')
            filestream.close()
            # run mg5 with the file generated
            if output is True:
                print('printing filename to screen')
            runcommand = 'cat ' + filename
            p = subprocess.run(runcommand, shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=mgloc + procloc)
            if output is True:
                print(p.stdout)
            runcommand = mgloc + procloc + '/bin/madevent ' + filename
            #print('runcommand=',runcommand)
            p = subprocess.Popen(runcommand, shell=True, text=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd=mgloc + procloc)
            for line in iter(p.stdout.readline, b''):
                pass
                if output is True:
                    print(line)
            if output is True:
                print(p.stdout)
                print(p.stderr)
            counter = counter + 1
    print('Done generating cross section')        
    return counter

def read_files(runnum, LambdasArray, m2, w2, m3, w3, process, nruns):
    if process in ProcLocation:
        procloc = ProcLocation[process]
    else:
        print('Process', process,'is not defined, exiting!')
    X = []
    Z = []
    XSEC = {}
    counter = 0
    for lams in LambdasArray:
        if counter > nruns:
            break
        #print(lams)
        lhe = 'run' + str(runnum) + '_m2_' + str(m2) + '_m3_' + str(m3) + '_w2_' + str(w2) + '_w3_' + str(w3) + '_' + '_'.join((lams)) + '/unweighted_events.lhe.gz'
        lhefile = MGLocation + procloc + 'Events/' + lhe
        #print('lhefile=', lhefile)
        #TestBool = True
        #if TestBool is False:
        if os.path.exists(lhefile) is False:
            print('Error, lhe file or summary file:', lhefile, 'does not exist!')
            exit()
        else:
            zgrepcommand = 'zgrep "Integrated weight" ' + lhefile
            p = subprocess.Popen(zgrepcommand, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, cwd='.')
            for line in iter(p.stdout.readline, b''):
                xsec = float(line.split()[5])
            #print(m2, m3, lams, xsec)
            #xsec = 0
            lams_tuple = []
            for mm in range(len(lams)):
                lams_tuple.append(float(lams[mm]))
            X.append(tuple(lams_tuple))
            Z.append(float(xsec))
            XSEC[tuple(lams_tuple)] = float(xsec)
            #print(X)
        counter = counter + 1
    #return np.transpose(X), Z, XSEC
    return X, Z, XSEC

def get_mg5_xsec(process, runnum, LambdasArray, k1, k2, k3, m2, w2, m3, w3, ecm=13):
    # check if file exists and run MG5 if not:
    drive_mg(process, runnum, MGLocation, round_sig(k1,4), round_sig(k2,4), round_sig(k3,4), [[str(round_sig(lam,4)) for lam in LambdasArray]], round_sig(m2,4), round_sig(w2,4), round_sig(m3,4), round_sig(w3,4), 1, 1,output=True, ecm=ecm)
    # read the xsec 
    X, Z, XSEC = read_files(runnum, [[str(round_sig(lam,4)) for lam in LambdasArray]], round_sig(m2,4), round_sig(w2,4), round_sig(m3,4), round_sig(w3,4), process, 1)
    return Z[0]


# TEST:
#Lambdas = [ 27.940871698413677 , 39.07032671865848 , 160.76112043143365 , -462.90567518636254 , -65.75926975945164 , 0.07061899381530284 , -0.19910092033686816 , -0.017181659533562252 , 160.76112043143365 ]
#k1 = 0.966
#k2 = 0.094
#k3 = 0.239
#m2 = 255
#m3 = 504
#w2 = 0.086
#w3 = 11
#xsec = get_mg5_xsec('hhh', RunNum, Lambdas, k1, k2, k3, m2, w2, m3, w3)
#print(xsec)
