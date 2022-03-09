#!/usr/bin/env python
"""
The Scripts is used to calculate GQS for single locus or given multiple locus.
"""
## import libraries
import argparse
import sys
import pandas as pd
import numpy as np
import os
import re
import subprocess
parser=argparse.ArgumentParser()
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from contextlib import suppress
import datetime

## fixed parameters
q1=0.68;q2=0.95;q3=0.99
nth=5
crit_cord=[(0.4,1),(1, 0.4)]
##

TODAY_YMD = datetime.datetime.today().strftime("Date: %Y-%m-%d \nTime: %H:%M:%S\n")
__version__ = '1.0.1'
DEFHEAD = "************************************************\n"
DEFHEAD += "* Genomic Quality Score (GQS)\n"
DEFHEAD += "* Version {V}\n".format(V=__version__)
DEFHEAD += "* Created by Swapnil Awasthi\n"
DEFHEAD += "***********************************************\n"
DEFHEAD += TODAY_YMD
DEFHEAD += "***********************************************\n"

parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--ifile', default=None, type=str,  help="Input filename (single/sumstats)")
parser.add_argument('--regs',  default=None, type=str,  help="Defined regions (with these columns in this order CHR,START and END)")
parser.add_argument('--r2_th', default=0.0,  type=float,help="SNPs below this LD will be excluded from the analysis (Default:0.0)")
parser.add_argument('--chrm',  default=None, type=str,  help="Chromosome number")
parser.add_argument('--chrm_h',default=None, type=str,  help="Column name for the chromosome in the input file")
parser.add_argument('--pval_h',default=None, type=str,  help="Column name for the Pvalue in the input file")
parser.add_argument('--snp_h', default=None, type=str,  help="Column name for the SNP identifier in the input file")
parser.add_argument('--pos_h', default=None, type=str,  help="Column name for the SNP position in the input file")
parser.add_argument('--ld_h',  default=None, type=str,  help="Column name for the LD in the input file")
parser.add_argument('--refG',  default=None, type=str,  choices=['genome1000-EUR','genome1000-EAS','HRC-EUR','HRC-EAS'], help='Reference Panel')
parser.add_argument('--addout',default=None, type=str,  help='prefix for the output files')

def sort_args(opts):
    """ Maps the column names"""
    cols_opts = [
        [opts['chrm_h'], 'CHR'],
        [opts['snp_h'],  'SNP'],
        [opts['pval_h'],'PVAL'],
        [opts['ld_h'],  'RSQR'],
        [opts['pos_h'],  'POS'],
    ]
    col_dict = {x[0]: x[1] for x in cols_opts if x[0] is not None}
    col_list=list(col_dict.keys())
    return [col_dict, col_list]

def read_ifile(ifile,col_list, col_dict):
    """Try to read the file correctly"""
    try:
        data = pd.read_csv(ifile, delim_whitespace=True, usecols=col_list)
        data=data.rename(columns=col_dict)
        return (data)
    except ValueError:
        raise ValueError('One of these columns did not matched to the header\n'+'\n'.join(col_list))

def verify_chr(data,chrm):
    """Interpret the chromosome number from the file or verify the entered one"""
    chr_list=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22']
    if chrm is not None:
        if chrm in chr_list:
            chrm=chrm
        else:
            raise ValueError('Not a valid chromosome number')

    else:
        nchrm=data.CHR.value_counts().rename_axis('unique_values').reset_index(name='counts')
        if nchrm.shape[0]>1:
            raise ValueError('Column contains more than one possible chromosome?')

        else:
            chrm=np.str(nchrm['unique_values'].values[0])
            if chrm in chr_list:
                chrm=chrm
            else:
                raise ValueError('Not a valid chromosome number')
    return chrm

class Logger(object):
    """Generate log files"""
    def __init__(self, fh):
        self.log_fh = open(fh, 'a')
    def log(self, msg):
        self.log_fh.write(msg)

def ld_filter(data,r2_th):
    """Filter in LD threshold return excluded and not included data"""
    P1_index = data[data['RSQR'] > r2_th].index
    P2_index = data[data['RSQR'] <=r2_th].index
    data_inc = data.loc[P1_index]
    data_exc = data.loc[P2_index]
    data_exc.to_csv(addout+'_exc.r2th.txt', index=None, sep='\t')
    data_inc.to_csv(addout+'_for.plot.txt', index=None, sep='\t')
    return data_inc, data_exc

def get_index(bol, data):
    """Try to find index SNP"""
    if bol is None:
        index=data.nsmallest(1,'PVAL', keep='first')
        return [index.SNP.values[0], index.PVAL.values[0]]
    else:
        index=data[data['RSQR']==1]
        index=index.nsmallest(1,'PVAL', keep='first')
        return [index.SNP.values[0], index.PVAL.values[0],index.PVAL_L.values[0]]

def calc_slope(cord):
    """Calculates slope of line"""
    x1=cord[0][0]; y1=cord[0][1]; x2=cord[1][0]; y2=cord[1][1]
    m=(y2-y1)/(x2-x1); i=y1-(m*x1)
    return [m,i]

def shortest_dist(line,p):
    """Calcualtes shortest distance between two points on 2D space"""
    m1=line[0]; i1=line[1]; x2=p[0]; y2=p[1]
    m2=-1/line[0]; i2=y2-(m2*x2)
    x_meet=(i2-i1)/(m1-m2); y_meet=m2*x_meet+i2
    dist=np.round(np.sqrt((x_meet-x2)**2+(y_meet-y2)**2),3)
    return dist

def color_code(index_snp, data):
    """Assigns color for the scatter plot"""
    # different quantile of absolute residuals
    P68 = data['2*RESID'].quantile(q1); P95 = data['2*RESID'].quantile(q2); P99 = data['2*RESID'].quantile(q3)
    # Colours codes 68%,95%,99%
    data['Colors'] = np.where(data['2*RESID'] >= P68, 'red', 'black')
    data.loc[(data['2*RESID'] >= P68) & (data['2*RESID'] < P95), 'Colors'] = 'green'
    data.loc[(data['2*RESID'] >= P95) & (data['2*RESID'] < P99), 'Colors'] = 'orange'
    data.loc[data['SNP'] == index_snp[0], 'Colors'] = 'magenta'  # index snp
    # directions
    data['SIGN'] = np.where(data.RESID > 0, '+', '-')
    data.loc[data.RESID == 0, 'SIGN'] = '0'
    ##areas
    data = data.sort_values(['PVAL_L_PRED'])
    data['Up68'] = np.round(data['PVAL_L_PRED'] + P68, 3); data['Lo68'] = np.round(data['PVAL_L_PRED'] - P68, 3)
    data['Up95'] = np.round(data['PVAL_L_PRED'] + P95, 3); data['Lo95'] = np.round(data['PVAL_L_PRED'] - P95, 3)
    data['Up99'] = np.round(data['PVAL_L_PRED'] + P99, 3); data['Lo99'] = np.round(data['PVAL_L_PRED'] - P99, 3)
    return data

def prep_gqs(data):
    """Estimates necessary parameters for calculating GQS"""
    data['PVAL_L'] = np.round(-np.log10(data['PVAL']), 3)       # negative log of p-values
    index_snp=get_index(True, data) # index snp
    b_cord=[(0, 0), (1, index_snp[2])]  # line of the orgin and index snp
    b_m, b_i = calc_slope(b_cord)

    data['PVAL_L_PRED'] = np.round(data['RSQR']*b_m + b_i, 3)   # expected 'PVAL_L' according to the above line
    data['RESID']   = np.round(data['PVAL_L'] - data['PVAL_L_PRED'], 3) # residuals
    data['2*RESID'] = np.round(np.sqrt(data['RESID']*data['RESID']), 3) # absolute values of residuals
    data=color_code(index_snp, data)

    # residuals/predicted p-value (i.e amount of signal lost or gain)
    data['M_RESID'] = np.round(data['RESID'] / data['PVAL_L_PRED'], 3)

    return index_snp, data

def data_sep(ndata):
    """Some necessary modification for the plots"""
    ndata_pos = ndata[(ndata.SIGN == '+')].shape[0]
    ndata_neg = ndata[(ndata.SIGN == '-')].shape[0] + ndata[(ndata.SIGN == '0')].shape[0]

    ndata68_pos = ndata[(ndata.SIGN == '+') & (ndata.Colors == 'black')].shape[0]
    ndata68_neg = ndata[(ndata.SIGN == '-') & (ndata.Colors == 'black')].shape[0] + ndata[(ndata.SIGN == '0')].shape[0]

    ndata95_pos = ndata[(ndata.SIGN == '+') & (ndata.Colors == 'green')].shape[0]
    ndata95_neg = ndata[(ndata.SIGN == '-') & (ndata.Colors == 'green')].shape[0]

    ndata99_pos = ndata[(ndata.SIGN == '+') & (ndata.Colors == 'orange')].shape[0]
    ndata99_neg = ndata[(ndata.SIGN == '-') & (ndata.Colors == 'orange')].shape[0]

    ndata10_pos = ndata[(ndata.SIGN == '+') & (ndata.Colors == 'red')].shape[0]
    ndata10_neg = ndata[(ndata.SIGN == '-') & (ndata.Colors == 'red')].shape[0]

    return ({'ndata_pos':ndata_pos,'ndata_neg':ndata_neg,'ndata68_pos':ndata68_pos,'ndata68_neg':ndata68_neg,\
             'ndata95_pos':ndata95_pos,'ndata95_neg':ndata95_neg,'ndata99_pos':ndata99_pos,\
             'ndata99_neg':ndata99_neg,'ndata10_pos':ndata10_pos, 'ndata10_neg':ndata10_neg})

def plot_ldvsPval(ndata):
    """Function to plot Scatter plot"""
    index = get_index(True, ndata)
    f, ax1 = plt.subplots(figsize=(10, 8))
    ax1.scatter(ndata['RSQR'],ndata['PVAL_L'], s=25, color=ndata['Colors'], label=None)
    ax1.scatter(1,index[2],s=40,marker='d',color='magenta',label=None)
    ax1.plot([1,0],[index[2],0],color='black')

    ##filling and colored line
    #green
    ax1.fill_between(ndata.RSQR, ndata.Lo68, ndata.Up68, color='#888888', alpha=0.30)
    lly68 = ndata.nsmallest(1, 'Lo68')['Lo68'].values; luy68 = ndata.nlargest(1, 'Lo68')['Lo68'].values
    llx68 = ndata.nsmallest(1, 'Lo68')['RSQR'].values; lux68 = ndata.nlargest(1, 'Lo68')['RSQR'].values
    ax1.plot([lux68, llx68], [luy68, lly68], color='green', linewidth=1)
    uly68 = ndata.nsmallest(1, 'Up68')['Up68'].values; uuy68 = ndata.nlargest(1, 'Up68')['Up68'].values
    ulx68 = ndata.nsmallest(1, 'Up68')['RSQR'].values; uux68 = ndata.nlargest(1, 'Up68')['RSQR'].values
    ax1.plot([uux68, ulx68], [uuy68, uly68], color='green', linewidth=1)
    #orange
    ax1.fill_between(ndata.RSQR, ndata.Lo95, ndata.Up95, color='#888888', alpha=0.20)
    lly95 = ndata.nsmallest(1, 'Lo95')['Lo95'].values; luy95 = ndata.nlargest(1, 'Lo95')['Lo95'].values
    llx95 = ndata.nsmallest(1, 'Lo95')['RSQR'].values; lux95 = ndata.nlargest(1, 'Lo95')['RSQR'].values
    ax1.plot([lux95, llx95], [luy95, lly95], color='orange', linewidth=1)
    uly95 = ndata.nsmallest(1, 'Up95')['Up95'].values; uuy95 = ndata.nlargest(1, 'Up95')['Up95'].values
    ulx95 = ndata.nsmallest(1, 'Up95')['RSQR'].values; uux95 = ndata.nlargest(1, 'Up95')['RSQR'].values
    ax1.plot([uux95, ulx95], [uuy95, uly95], color='orange', linewidth=1)
    #red
    ax1.fill_between(ndata.RSQR, ndata.Lo99, ndata.Up99, color='#888888', alpha=0.10)
    lly99 = ndata.nsmallest(1, 'Lo99')['Lo99'].values; luy99 = ndata.nlargest(1, 'Lo99')['Lo99'].values
    llx99 = ndata.nsmallest(1, 'Lo99')['RSQR'].values; lux99 = ndata.nlargest(1, 'Lo99')['RSQR'].values
    ax1.plot([lux99, llx99], [luy99, lly99], color='red', linewidth=1)
    uly99 = ndata.nsmallest(1, 'Up99')['Up99'].values; uuy99 = ndata.nlargest(1, 'Up99')['Up99'].values
    ulx99 = ndata.nsmallest(1, 'Up99')['RSQR'].values; uux99 = ndata.nlargest(1, 'Up99')['RSQR'].values
    ax1.plot([uux99, ulx99], [uuy99, uly99], color='red', linewidth=1)
    sep=data_sep(ndata)

    columns = ('Interval', '-SNP', '+SNP')
    cell_text = [['0',     str(sep['ndata_neg']),   str(sep['ndata_pos'])], ['0-68', str(sep['ndata68_neg']), str(sep['ndata68_pos'])],
                 ['68-95', str(sep['ndata95_neg']), str(sep['ndata95_pos'])], ['95-99', str(sep['ndata99_neg']), str(sep['ndata99_pos'])],
                 ['99-100',str(sep['ndata10_neg']), str(sep['ndata10_pos'])]]

    the_table2 = plt.table(cellText=cell_text, colLabels=columns, loc='upper left', cellLoc='center',
                           colWidths=[0.07, 0.057, 0.060], \
                           colColours=['whitesmoke'] * 3,
                           cellColours=[['whitesmoke'] * 3, ['grey'] * 3, ['green'] * 3, ['orange'] * 3, ['red'] * 3])

    the_table2.auto_set_font_size(False)
    the_table2.set_fontsize(11)
    the_table2.scale(1.2, 1.2)
    plt.suptitle('Index SNP: ' + index[0], fontsize=20)
    plt.suptitle('Index SNP: ' + index[0], fontsize=20)
    plt.title('r2 > '+str(args.r2_th), fontsize=18)

    plt.ylim(0, index[2] + 1)
    plt.xlim(0, 1.01)

    plt.xlabel("r^2", fontsize=14)
    plt.ylabel("Observed (-logP)", fontsize=14)
    plt.savefig(addout+'_R2vsPl.plot.pdf', dpi=400)
    plt.savefig(addout+'_R2vsPl.plot.png', dpi=400)

def plot_hist(ndata):
    """Function to plot histogram"""
    n = ndata.shape[0]
    f, ax2 = plt.subplots(figsize=(10, 8))
    P68 = ndata['2*RESID'].quantile(0.68); P95 = ndata['2*RESID'].quantile(0.95)
    P99 = ndata['2*RESID'].quantile(0.99)
    index = get_index(True, ndata)
    ax2.hist(ndata['RESID'], bins=150, color='grey')
    ax2.axvline(x=0, color='black');    ax2.axvline(x=P68, color='green')
    ax2.axvline(x=-P68, color='green'); ax2.axvline(x=P95, color='orange')
    ax2.axvline(x=-P95, color='orange');ax2.axvline(x=P99, color='red')
    ax2.axvline(x=-P99, color='red')
    sep=data_sep(ndata)
    L1_per_pos = str(np.round((sep['ndata_pos'] / n) * 100, 1))
    L1_per_neg = str(np.round((sep['ndata_neg'] / n) * 100, 1))
    L168_per_pos = str(np.round((sep['ndata68_pos'] / n) * 100, 1))
    L168_per_neg = str(np.round((sep['ndata68_neg'] / n) * 100, 1))
    L195_per_pos = str(np.round((sep['ndata95_pos'] / n) * 100, 1))
    L195_per_neg = str(np.round((sep['ndata95_neg'] / n) * 100, 1))
    L199_per_pos = str(np.round((sep['ndata99_pos'] / n) * 100, 1))
    L199_per_neg = str(np.round((sep['ndata99_neg'] / n) * 100, 1))

    L1100_per_pos = str(np.round((sep['ndata10_pos'] / n) * 100, 1))
    L1100_per_neg = str(np.round((sep['ndata10_neg'] / n) * 100, 1))

    columns = ('Interval', '-SNP (%)', '+SNP (%)')

    cell_text = [['0', L1_per_neg, L1_per_pos], ['0-68', L168_per_neg, L168_per_pos],
                 ['68-95', L195_per_neg, L195_per_pos], \
                 ['95-99', L199_per_neg, L199_per_pos], ['99-100', L1100_per_neg, L1100_per_pos]]

    the_table2 = plt.table(cellText=cell_text, colLabels=columns, loc='upper left', cellLoc='center',
                           colWidths=[0.075, 0.09, 0.09], \
                           colColours=['whitesmoke'] * 3,
                           cellColours=[['whitesmoke'] * 3, ['grey'] * 3, ['green'] * 3, ['orange'] * 3, ['red'] * 3])

    the_table2.auto_set_font_size(False)
    the_table2.set_fontsize(12)
    the_table2.scale(1.2, 1.2)
    plt.suptitle('Index SNP: ' + index[0], fontsize=20)
    plt.xlabel("RESID", fontsize=14)
    plt.ylabel("Number of SNPs", fontsize=14)
    plt.savefig(addout+'_Mresid.hist.pdf', dpi=400)
    plt.savefig(addout+'_Mresid.hist.png', dpi=400)

def mod_gqs(data):
    """Modification for GQS plot"""
    data_index = data[data.M_RESID <= 0.0].index
    ndata = data.loc[data_index]
    #This is critcial line
    c1_m,c1_i =calc_slope(crit_cord)
    line = (c1_m, c1_i)
    pw = (1.0, 1.0)  # cordinates for the worst case
    dist = shortest_dist(line, pw)
    ndata['dist_critical'] = ndata.apply(lambda row: shortest_dist(line, (row['RSQR'], -row['M_RESID'])), axis=1)
    ndata['gGQS'] = np.round((dist - ndata['dist_critical']) / dist, 3)

    ndata['M_RESID_pred'] = np.round((c1_m * ndata['RSQR'] + c1_i), 3)
    ndata['Rdist'] = np.round(abs(ndata['M_RESID']) - ndata['M_RESID_pred'], 3)
    ndata['Colors2'] = np.where(ndata['Rdist'] > 0.0, 'red', 'black')
    ndata.loc[(ndata['Rdist'] == 0.00), 'Colors2'] = 'yellow'

    return ndata


def plot_gqs(data):
    "Function to plot GQS"
    ndata=mod_gqs(data)
    f, ax2 = plt.subplots(figsize=(10, 8))
    ax2.scatter(ndata['RSQR'], -ndata['M_RESID'], s=30, edgecolor='black', color=ndata['Colors2'])
    plt.grid(True)

    nrdata_index = ndata[ndata.Colors2 == "red"].index
    nrdata = ndata.loc[nrdata_index]

    if ndata[ndata.RSQR>0.40].shape[0]<=1:
        gqs = -1
    elif nrdata.shape[0] == 0:
        gqs = 1
    else:

        gqs = nrdata.gGQS.min()
        gqs_index = nrdata[nrdata.gGQS == gqs].index
        gqs_df= nrdata.loc[gqs_index]
        lgqs = gqs_df.SNP.values[0]

        c2_i = np.abs(gqs_df['M_RESID'].values[0]) - 1 * gqs_df['RSQR'].values[0]
        x_corr = (c2_i - 1.4) /-2
        y_corr = -1.0 * x_corr + 1.4
        ax2.plot([gqs_df['RSQR'].values[0], x_corr], [np.abs(gqs_df['M_RESID'].values[0]), y_corr], linestyle='--',
                 color='coral')
        ax2.text(gqs_df['RSQR'].values[0], np.abs(gqs_df['M_RESID'].values[0]), lgqs)
        gqs = np.round(gqs, 3)
    ax2.plot([0.4, 1], [1, 0.4], linestyle='-', color='black')
    ax2.plot([1, 0.7], [1, 0.7], linestyle='--', color='black')

    plt.ylim(0, 1.01)
    plt.xlim(0, 1.01)
    plt.text(.90, .95, str(nrdata.shape[0]), fontsize=20, color='b')

    if (gqs == -1):
        plt.text(.40, 1.1, 'GQS: ' + str(gqs), fontsize=20, color='black', bbox=dict(boxstyle="round", fc='red'))
    elif(gqs == 1):
        plt.text(.40, 1.1, 'GQS: ' + str(gqs), fontsize=20, color='black', bbox=dict(boxstyle="round", fc='green'))
    else:
        plt.text(.40, 1.1, 'GQS: ' + str(gqs), fontsize=20, color='black', bbox=dict(boxstyle="round", fc='red'))

    plt.title('r^2 > '+str(args.r2_th), fontsize=18)
    plt.xlabel("r^2", fontsize=14)
    plt.ylabel("signal lost", fontsize=14)

    data.to_csv(addout+'_full.plot.txt',index=None, sep='\t')
    ndata.to_csv(addout + '_gqs.plot.txt', index=None, sep='\t')
    plt.savefig(addout+'_GqsVis.plot.pdf', dpi=400)
    plt.savefig(addout+'_GqsVis.plot.png', dpi=400)
    return [str(gqs),str(nrdata.shape[0])]

def cal_ld(ref_path, plink_path, chrm, data):
    """Calculate LD for each SNP w.r.t index using PLINK: Might need some work; Currently look for patter 'chr1'"""
    content=os.listdir(ref_path)
    r1 = re.compile(".*chr"+chrm+"[.|_].*fam") ##search patern--> chr1#
    fam_file = list(filter(r1.match, content))
    pref='.'.join(fam_file[0].split('.')[:-1])
    index=get_index(None, data)
    cmd1=plink_path+'plink --bfile '+ref_path+pref+' --r2 --ld-snp '+index[0]+' --ld-window-kb 1000 --ld-window 99999 --ld-window-r2 0 --out '+addout
    proc1=subprocess.Popen(cmd1, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    (out, err) = proc1.communicate()
    if not err.decode("utf-8").strip():
        pass
    else:
       raise ValueError('Are you sure the chromosome number is correct for this region/file')

    ld_info=pd.read_table(addout+'.ld', delim_whitespace=True)
    data_merge=pd.merge(data,ld_info,left_on=['SNP'],right_on=['SNP_B'],how='left')
    data_merge=data_merge[['SNP','PVAL','R2']]; data_merge.columns=['SNP','PVAL','RSQR']
    data_excl=data_merge[data_merge['RSQR'].isnull()]
    data_plot=data_merge[data_merge['RSQR'].notnull()]
    data_excl.to_csv(addout+'_exc.mis.ref.txt', index=None, sep='\t')
    return data_plot, data_excl

def create_config(name):
    """Setup configuration file: Might need some work for on other system"""
    cwd = os.getcwd()
    if not os.path.exists(cwd+'/'+name):
        open(cwd+'/'+name,'w+')
        #os.mknod(cwd+'/'+name)
        plink_loc = input("Enter the location of plink tool: ")
        genome1000EUR_loc = input("Enter the location of genome1000-EUR: ")
        genome1000EAS_loc = input("Enter the location of genome1000-EAS: ")
        HRCEUR_loc = input("Enter the location of HRC-EUR: ")
        HRCEAS_loc = input("Enter the location of HRC-EAS: ")
        lw1='plink_loc '+plink_loc+'\ngenome1000-EUR_loc '+genome1000EUR_loc+'\ngenome1000-EAS_loc '+\
            genome1000EAS_loc+'\nHRC-EUR_loc '+HRCEUR_loc+'\nHRC-EAS_loc '+HRCEAS_loc
        fw1=open(name,'a');fw1.write(lw1);fw1.close()
        sys.exit('Restart the script!')
    else:
        gqs_config={}
        with open(name, 'r') as f:
              for line in f:
                  line=line.strip().split(' ')
                  gqs_config[line[0]]=line[1]
        f.close()
        return gqs_config

def read_locus(data, chrm, start, end):
    """Seperate defined locus from summary statistic"""
    locus_data=data[(data['CHR']==int(chrm)) & (data['POS']>=int(start)) & (data['POS']<=int(end))]
    if locus_data.shape[0]==0:
        None
    else:
        return locus_data
def run_steps(data):
    """Sub function that call plot functions"""
    data_inc, data_exc2 = ld_filter(data, opts['r2_th'])
    if data_inc.shape[0] < nth:
        raise ValueError('The number of SNPs are very few after merging with the reference.\nPlease '
                         'match the SNP identifier with the reference')
    index, ndata = prep_gqs(data_inc)
    plot_ldvsPval(ndata)
    plot_hist(ndata)
    gqs = plot_gqs(ndata)
    out_res1='File_Name\tCHR\tIndex_SNP\tPvalue\tGQS\tOutliers\n'
    out_res2 = args.ifile + '\t' + chrm + '\t' + index[0] + '\t' + str('{:.3e}'.format(index[1])) + '\t' + gqs[
        0] + '\t' + gqs[1] + '\n'
    log_ob.log(out_res1)
    log_ob.log(out_res2)
    print(out_res1)
    print(out_res2)

def run_LD(data,chrm):
    stat = create_config('gqs_config')
    data, data_exc1 = cal_ld(stat[args.refG + '_loc'] + '/', stat['plink_loc'] + '/', chrm, data)
    if data.shape[0] < nth:
        raise ValueError('The number of SNPs are very few after merging with the reference.\nPlease '
                         'match the SNP identifier with the reference')
    return data, data_exc1
if __name__ == '__main__':
    args=parser.parse_args()
    if args.ifile is None:
        raise ValueError('The --ifile flag is required.')
    elif args.snp_h is None or args.pval_h is None or args.addout is None:
        raise ValueError('--snp_h, --pval_h and --addout flags are required.')
    elif args.chrm_h is None and args.chrm is None:
        raise ValueError('Either --chrm_h and --chrm flag is required.')
    elif args.ld_h is None and args.refG is None:
        raise ValueError('Either --ld_h and --refG flag is required.')
    if args.regs is not None:
        if args.pos_h is None:
            raise ValueError('--pos_h position column required to run on full sumstats')


    addout = args.addout
    opts = vars(args)
    [col_dict, col_list] = sort_args(opts)

    if args.regs is None:
        data = read_ifile(opts['ifile'], col_list, col_dict)
        chrm = verify_chr(data, args.chrm)
        if args.ld_h is None:
            data, data_exc1=run_LD(data,chrm)
        log_ob = Logger(addout + '.log')
        log_ob.log(DEFHEAD)
        run_steps(data)
    else:
        data = read_ifile(opts['ifile'], col_list, col_dict)
        with open(args.regs, 'r') as f1:
            addout_cont=addout; args.ifile_cont=args.ifile
            for line1 in f1:
                line1 = line1.strip().split()
                chrm = line1[0]; start = line1[1]; end = line1[2]
                addout = addout_cont + "_" + chrm + "_" + start + "_" + end
                args.ifile=args.ifile_cont + "_" + chrm + "_" + start + "_" + end
                locus_data=read_locus(data, chrm, start, end)

                log_ob = Logger(addout + '.log')
                log_ob.log(DEFHEAD)
                if locus_data is None:
                    out_res="No SNPs in that range"
                    log_ob.log(out_res)
                    print("No SNPs in that range")
                else:
                    locus_data, data_exc1 = run_LD(locus_data)
                    run_steps(locus_data)


