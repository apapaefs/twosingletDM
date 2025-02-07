import math
import sys
import gzip

# this python code generates an input file for MG5 to run and get the LO cross section for a specific process,
# while varying the mass.
# in this case, the chosen masses are read from the first column of an input file
# that corresponds to the Yellow Report cross section

#  generate the file:
generate = False

# read the results:
readresults = True

# file to read the masses from:
massfile = 'higgsXS_YR4_13TeV_NNLONNLL.txt'

# output script for madgraph
outputfile = 'gg_iota0_massvar.script'
# output file with cross sections
outputfileXS = 'higgsXS_136TeV_LO.txt'
# the subdirectory of the desired process (include slash at end)
procdir = 'gg_iota0/'
# name of the mass variable to vary
massvariable = 'Miota'

# madgraph directory to read results from and write xsec file (in pb) (include slash at end)
mgdir = '/mnt/ssd2/Projects/TwoSinglet/MG5_aMC_v2_9_15/'

# read the masses from massfile and put them in an array:
if generate is True:
    massarray = []
    f = open(massfile, "r")
    for line in f:
        massarray.append(float(line.split()[0]))
    print('Masses found in', massfile, '=', massarray)
    f.close()
    # generate the script:
    fout = open(outputfile, "w")
    fout.write('launch ' + str(procdir) + ' -i\n\n')
    for mass in massarray:
        fout.write('launch run_m' + str(mass) + ' --accuracy=0.25 --points=300 --iterations=1\n')
        fout.write('0\n')
        fout.write('set ' + str(massvariable) + ' ' + str(mass) + '\n')
        fout.write('0\n')

if readresults is True:
    massarray = []
    f = open(massfile, "r")
    for line in f:
        massarray.append(float(line.split()[0]))
    print('Masses found in', massfile, '=', massarray)
    f.close()
    fxs = open(outputfileXS, "w")
    for mass in massarray:
        resultsfile = mgdir + procdir + 'Events/' + 'run_m' + str(mass) + '/run_m' + str(mass) + '_tag_1_banner.txt'
        fin = open(resultsfile, "r")
        for line in fin:
            if 'Integrated weight' in line:
                xs = line.split()[-1]
                fxs.write(str(mass) + '\t' + str(xs) + '\n')
    fxs.close()
