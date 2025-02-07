import random
import math
import os
from sys import argv
from generate_trsm_info import * # TRSM info generator (branching ratios, mixing matrices, etc.)
from trsm_kstoalphas import * # Convert from k1, k2, k3 to a12, a13, a23
from prettytable import PrettyTable
from datetime import date
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker
from matplotlib.ticker import MultipleLocator
from tqdm import tqdm
from matplotlib import ticker, cm
from matplotlib.lines import Line2D
import scipy
from scipy import stats
from tqdm import tqdm
####################

# Energy for xsec calculation
Energy = 13.6

# SM XSEC at Energy [pb]
xsec_sm = {}
xsec_sm[13] = 3.641e-05
xsec_sm[13.6] = 4.029e-05

# TAG for RUN to plot:
#ini_seed = 1234
#DATE = 20230623
#DATE = 20230702
DATE = 20230825
ini_seed = 9999

# Directory for the output files to be read
OutputDir = 'output/'

# directory for plots:
plotdir = 'plots/'

# filename to be processed
# my files, automatically generated
#RunTag = str(Energy) + '-' + str(DATE) + '-' + str(ini_seed) + '-True'
#filename = OutputDir + 'trsm_points_' + RunTag + '.dat'
# Gilberto's scan:
RunTag = 'GTX_test'
filename = OutputDir + 'All_files_final.txt'

####################
# functions
####################

# choose the next colour -- for plotting
colors = [ 'green', 'orange', 'red', 'magenta', 'blue', 'cyan', 'black', 'brown', 'violet'] # 9 colours

ccount = 0
def next_color():
    global ccount
    color_chosen = colors[ccount]
    if ccount < 8:
        ccount = ccount + 1
    else:
        ccount = 0    
    return color_chosen

# do not increment colour in this case:
def same_color():
    global ccount
    color_chosen = colors[ccount-1]
    return color_chosen

# reset the colour counter
def reset_color():
    global ccount
    ccount = 0

# change the colour scheme   
def change_color_scheme():
    global colors
    colors = ['violet', 'blue', 'green', 'orange', 'red', 'black', 'magenta', 'brown', 'cyan']

# change the colour scheme   
def change_color_scheme2():
    global colors
    colors = ['brown', 'cyan', 'green', 'orange', 'red', 'black', 'magenta', 'brown', 'cyan']

# print the parameter point info info:
def print_info(vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3, fxsec136):
    
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
    tbl.add_row(['fxsec136', fxsec136])
    print(tbl)
    #print('\n')


energyvar = ['m_2','m_3', 'v_x', 'v_s']



# plot viable points on the m2-m3 plane:
def scatter_plot(data, var1, var2):
    plot_type = 'scatter_plot_' + str(var1).replace('{','').replace('}','').replace('_','') + '_' + str(var2).replace('{','').replace('}','').replace('_','') + '_' + RunTag.replace('-True','')
    # the following labels are in LaTeX, but instead of a single slash, two "\\" are required.
    xlab = var1 # the x label
    ylab = var2 # the y label
    # set up matplotlib
    gs = gridspec.GridSpec(4, 4)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid(False)
    # create plot title 
    plot_title = 'Viable points (RGEs, Other Theory, Expt. Constraints)'
    ax.set_title(plot_title)

    # set the colors for the points
    colors = []
    alphas = []
    for d in data['fxsec']:
        if d > 10. and d < 20.:
            colors.append('blue')
            alphas.append(0.5)
        elif d > 20. and d < 50.:
            colors.append('green')
            alphas.append(0.5)
        elif d > 50. and d < 100.:
            colors.append('red')
            alphas.append(0.5)
        elif d > 100.:
            colors.append('magenta')
            alphas.append(0.5)
        else:
            colors.append('black')
            alphas.append(0.2)
    # and the custom legend corresponding to these:
    legend_elements  = [Line2D([0], [0], marker='o', color='w', label=r'$< \times 10$', markerfacecolor='black', markersize=8, alpha=0.8), Line2D([0], [0], marker='o', color='w', label=r'$\times (10-20)$', markerfacecolor='blue', markersize=8, alpha=0.8), Line2D([0], [0], marker='o', color='w', label=r'$\times (20-50)$', markerfacecolor='green', markersize=8, alpha=0.8), Line2D([0], [0], marker='o', color='w', label=r'$ \times (50-100)$', markerfacecolor='red', markersize=8, alpha=0.8), Line2D([0], [0], marker='o', color='w', label=r'$> \times 100$', markerfacecolor='magenta', markersize=8, alpha=0.8)]
    
    ax.legend(handles=legend_elements, loc='upper right',framealpha=0.6, fancybox=True, title=r'$\times \sigma_\mathrm{SM}(gg\rightarrow hhh)@' + str(Energy) + '$ TeV')
    # scatter plot
    plt.scatter(data[var1], data[var2], s=data['fxsec'], c=colors, alpha=alphas)

    # set the x and y labels
    if var1 in energyvar:
        xlab = '$' + var1 + '$ [GeV]'
    else:
        xlab = '$' + var1 + '$'
    if var2 in energyvar:
        ylab = '$' + var2 + '$ [GeV]'
    else:
        ylab = '$' + var2 + '$'
    ax.set_ylabel(ylab, fontsize=20)
    ax.set_xlabel(xlab, fontsize=20)
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    # write out the file
    infile = plot_type + '.dat'
    print('---')
    print('output in', plotdir + infile.replace('.dat','.pdf'))
    plt.savefig(plotdir + infile.replace('.dat','.pdf'), bbox_inches='tight')
    plt.close(fig)


# plot viable points on a chosen plane, no marker size according to xsec
def scatter_plot_simple(data, var1, var2):
    plot_type = 'scatter_plot_simple_' + str(var1).replace('{','').replace('}','').replace('_','') + '_' + str(var2).replace('{','').replace('}','').replace('_','') + '_' + RunTag.replace('-True','')
    # the following labels are in LaTeX, but instead of a single slash, two "\\" are required.
    xlab = var1 # the x label
    ylab = var2 # the y label
    # set up matplotlib
    gs = gridspec.GridSpec(4, 4)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid(False)
    # create plot title 
    plot_title = 'Viable points (RGEs, Other Theory, Expt. Constraints)'
    ax.set_title(plot_title)

    #ax.legend(loc='upper right',framealpha=0.6, fancybox=True, title=r'$\times \sigma_\mathrm{SM}(gg\rightarrow hhh)@' + str(Energy) + '$ TeV')
    # scatter plot
    plt.scatter(data[var1], data[var2], s=15, c='red', alpha=0.5)

    # set the x and y labels
    if var1 in energyvar:
        xlab = '$' + var1 + '$ [GeV]'
    else:
        xlab = '$' + var1 + '$'
    if var2 in energyvar:
        ylab = '$' + var2 + '$ [GeV]'
    else:
        ylab = '$' + var2 + '$'

    # special cases here:
    if var2 == 'fxsec':
        ylab = r'$\times \sigma_\mathrm{SM}(gg\rightarrow hhh)@' + str(Energy) + '$ TeV'
    if var1 == 'resxsec':
        xlab = r'fraction from resonant: $pp \rightarrow h_3 \rightarrow h_2 h_1 \rightarrow h_1 h_1 h_1$'
    if var1 == 'k123xk112sq':
        xlab = r'$\sqrt{|\lambda_{123}\times \lambda_{112}|}$ [GeV]'

    # write out the file
    # log scale?
    ylog = True # whether to plot y in log scale
    xlog = False # whether to plot x in log scale
    # choose x and y log scales
    if ylog:
        ax.set_yscale('log')
    else:
        ax.set_yscale('linear')
    if xlog:
        ax.set_xscale('log')
    else:
        ax.set_xscale('linear')
        
    ax.set_ylabel(ylab, fontsize=15)
    ax.set_xlabel(xlab, fontsize=15)
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    # write out the file
    infile = plot_type + '.dat'
    print('---')
    print('output in', plotdir + infile.replace('.dat','.pdf'))
    plt.savefig(plotdir + infile.replace('.dat','.pdf'), bbox_inches='tight')
    plt.close(fig)


# plot viable points on a chosen plane, no marker size according to xsec
def density_plot_simple(data, var1, var2):
    plot_type = 'density_plot_simple_' + str(var1).replace('{','').replace('}','').replace('_','') + '_' + str(var2).replace('{','').replace('}','').replace('_','') + '_' + RunTag.replace('-True','')
    # the following labels are in LaTeX, but instead of a single slash, two "\\" are required.
    xlab = var1 # the x label
    ylab = var2 # the y label
    # set up matplotlib
    gs = gridspec.GridSpec(4, 4)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid(False)
    # create plot title 
    plot_title = 'Viable points (RGEs, Other Theory, Expt. Constraints)'
    ax.set_title(plot_title)

    # name the data as x and y
    x = data[var1]
    y = data[var2]
    
    # # Calculate the point density
    xy = np.vstack([x,y])
    z = scipy.stats.gaussian_kde(xy)(xy)

    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    
    #ax.legend(loc='upper right',framealpha=0.6, fancybox=True, title=r'$\times \sigma_\mathrm{SM}(gg\rightarrow hhh)@' + str(Energy) + '$ TeV')
    # scatter plot with density info:
    ax.scatter(x, y, c=z, s=50, edgecolor=['none'])
    #plt.scatter(data[var1], data[var2], s=15, c='red', alpha=0.5)

    # log scale?
    ylog = True # whether to plot y in log scale
    xlog = False # whether to plot x in log scale

    # set the x and y labels
    if var1 in energyvar:
        xlab = '$' + var1 + '$ [GeV]'
    else:
        xlab = '$' + var1 + '$'
    if var2 in energyvar:
        ylab = '$' + var2 + '$ [GeV]'
    else:
        ylab = '$' + var2 + '$'

    # special cases here:
    if var2 == 'fxsec':
        ylab = r'$\times \sigma_\mathrm{SM}(gg\rightarrow hhh)@' + str(Energy) + '$ TeV'
    if var1 == 'resxsec':
        xlab = r'fraction from resonant: $pp \rightarrow h_3 \rightarrow h_2 h_1 \rightarrow h_1 h_1 h_1$'
    if var1 == 'k123xk112sq':
        xlab = r'$\sqrt{|\lambda_{123}\times \lambda_{112}|}$ [GeV]'

    # write out the file

    # choose x and y log scales
    if ylog:
        ax.set_yscale('log')
    else:
        ax.set_yscale('linear')
    if xlog:
        ax.set_xscale('log')
    else:
        ax.set_xscale('linear')
        
    ax.set_ylabel(ylab, fontsize=15)
    ax.set_xlabel(xlab, fontsize=15)
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    # write out the file
    infile = plot_type + '.dat'
    print('---')
    print('output in', plotdir + infile.replace('.dat','.pdf'))
    plt.savefig(plotdir + infile.replace('.dat','.pdf'), bbox_inches='tight')
    plt.close(fig)    


def plot_ks(data):
    plot_type = 'scatter_plot_k_vs_M_' + RunTag.replace('-True','')
    # set up matplotlib
    gs = gridspec.GridSpec(4, 4)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid(False)
    # create plot title 
    plot_title = 'Viable points (RGEs, Other Theory, Expt. Constraints)'
    ax.set_title(plot_title)
    ylab='$k_i, i=2,3$ or $1-k_1$'
    xlab='$m_i$ [GeV]'
    ax.set_ylabel(ylab, fontsize=20)
    ax.set_xlabel(xlab, fontsize=20)
    # scatter plot
    for i in range(1,4):
        if i == 1:
            plt.scatter(data['m_' + str(i)],1-data['k_' + str(i)], s=data['fxsec'], c=next_color(), alpha=0.9, label='$1-k_1$')
        else:
            plt.scatter(data['m_' + str(i)],data['k_' + str(i)], s=data['fxsec'], c=next_color(), alpha=0.9,label='$k_' + str(i) + '$')
    reset_color()
    ax.legend(loc='lower right',framealpha=1.0, fancybox=True)

    # write out the file
    # log scale?
    ylog = True # whether to plot y in log scale
    xlog = False # whether to plot x in log scale
    # choose x and y log scales
    if ylog:
        ax.set_yscale('log')
    else:
        ax.set_yscale('linear')
    if xlog:
        ax.set_xscale('log')
    else:
        ax.set_xscale('linear')
    ax.xaxis.set_major_locator(ticker.AutoLocator())
    ax.xaxis.set_minor_locator(ticker.AutoMinorLocator())

    infile = plot_type + '.dat'
    
    print('---')
    print('output in', plotdir + infile.replace('.dat','.pdf'))
    plt.savefig(plotdir + infile.replace('.dat','.pdf'), bbox_inches='tight')
    plt.close(fig)

    
############################
# MAIN FUNCTION HERE
############################


# MAXPOINTS for testing:
#MAXPOINTS = 1000
    
# open file and read in data:
df = pd.read_table(filename, sep ='\t', header=None, names=['m_2', 'm_3', 'v_s', 'v_x', 'a_{12}', 'a_{13}','a_{23}', 'fxsec', 'resxsec'])#,nrows=MAXPOINTS)

#print(df)

# generate more information:
k1arr = []
k2arr = []
k3arr = []
m1arr = []
k123sqarr = []
k112sqarr = []
k123xk112sqarr = []
pbar = tqdm(total=df.size)
for data in df.itertuples():
    vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3, h1_BRs, h2_BRs, h3_BRs, xs136_lo_h1, xs136_lo_h2, xs136_lo_h3 = generate_lams(ini_seed, data.m_2, data.m_3, data.v_s, data.v_x, data._5, data._6, data._7, False)
    #if data.fxsec > 50.:
    #    print_info(vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3, data.fxsec)
    k1arr.append(k1)
    k2arr.append(k2)
    k3arr.append(k3)
    m1arr.append(125.)
    k123sqarr.append(K123**2)
    k112sqarr.append(K112**2)
    k123xk112sqarr.append(math.sqrt(abs(K123 * K112)))   
    pbar.update(1)
    

#for i, data in df.iterrows():
    #print(data['m_2'])
    #print(data['m_2'][i], data['m_3'][i], data['v_s'][i], data['v_x'][i], data['a_{12}'][i], data['a_{13}'][i], data['a_{23}'][i])
#    vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3, h1_BRs, h2_BRs, h3_BRs, xs136_lo_h1, xs136_lo_h2, xs136_lo_h3 = generate_lams(ini_seed, data['m_2'], data['m_3'], data['v_s'], data['v_x'], data['a_{12}'], data['a_{13}'], data['a_{23}'], False)
#    if df['fxsec'][i] > 50.:
#        print_info(vs, vx, M2, M3, a12, a13, a23, w1, w2, w3, K111, K112, K113, K123, K122, K1111, K1112, K1113, K133, k1, k2, k3, df['fxsec'][i])
#    k1arr.append(k1)
#    k2arr.append(k2)
#    k3arr.append(k3)
#    m1arr.append(125.)
#    k123sqarr.append(K123**2)
#    k112sqarr.append(K112**2)
#    k123xk112sqarr.append(math.sqrt(abs(K123 * K112)))
    
df['k_1'] = k1arr
df['k_2'] = k2arr
df['k_3'] = k3arr
df['m_1'] = m1arr
df['k123sq'] = k123sqarr
df['k112sq'] = k112sqarr
df['k123xk112sq'] = k123xk112sqarr

# print the table
#print(df)
#print(df['fxsec'])

dfsorted = df.sort_values(by=['fxsec'], ascending=False)
print(dfsorted[:12])

# plots:
scatter_plot(df, 'm_2', 'm_3')
scatter_plot_simple(df, 'resxsec', 'fxsec')
density_plot_simple(df, 'resxsec', 'fxsec')
density_plot_simple(df, 'k123xk112sq', 'fxsec')

plot_ks(df)


