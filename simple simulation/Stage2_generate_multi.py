import sys
import os
from numpy import random 
import random as rn
import math
from scipy.stats import expon
from collections import Counter
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


global CL
CL = 64444167
global MaMC
MaMC = 20
global MiMC
MiMC = 1
global ConstMulti_for_RNAmulti
ConstMulti_for_RNAmulti = 5
global ConstMulti_for_DNAmulti
ConstMulti_for_DNAmulti = 6
global ConstMulti_for_DNAmultiRNAmulti
ConstMulti_for_DNAmultiRNAmulti = 4

global proportion
proportion = "1:2:1"#rna_multi:dna_multi:rna_dna_multi

global exp_scale
lambd = -math.log(0.01)/MaMC
exp_scale = 1/lambd 

outputPath = sys.argv[1] #"/home/berdnikovich/try_mHiC/result_simulation_simmple/new_generation/sim1/s4"#sys.argv[1]#'/home/berdnikovich/try_mHiC/simulation/10_sim/sim1/s4
libName = '2stage{}'.format(MaMC)
picpath = sys.argv[2] #"/home/berdnikovich/try_mHiC/result_simulation_simmple/new_generation/sim1/s4"#sys.argv[2]

###used files
global logfile
logfile = os.path.join(outputPath, libName+".fithic.log")
global multitrue
multitrue = os.path.join(outputPath, 'Multi_true_coo')
global output1
output1 = os.path.join(outputPath, 'Multi_tf')

###
global pic1
pic1 = os.path.join(picpath, 'Stage2_Multi_false{}.png'.format(MaMC))
global pic2
pic2 = os.path.join(picpath, 'Stage2_Multi_false_true{}.png'.format(MaMC))
global pic3
pic3 = os.path.join(picpath,'Stage2_Multi_counts{}.png'.format(MaMC))

with open(logfile, 'w') as log:
        log.write("\n\nStart gene multi cont\n")
        log.write("\ndist_func_multi: uniform\n")
        log.write("\nNumber of pseudo pairs:\n")
        log.write("\nRNA multi{}\n".format(ConstMulti_for_RNAmulti))
        log.write("\nDNA multi{}\n".format(ConstMulti_for_DNAmulti))
        log.write("\nRNA multi - DNA multi{}\n".format(ConstMulti_for_DNAmultiRNAmulti))
        
        
with open(multitrue, 'r') as u_t_f:
    uni_true = []
    for ut in u_t_f:
        ut_ = ut.strip().split('\t')
        uni_true.append([int(ut_[0]), int(ut_[1])])
UtC = len(uni_true)


shares = [int(i) for i in proportion.split(":")]
rna_m_c = int((UtC/sum(shares))*shares[0])
dna_m_c = int((UtC/sum(shares))*shares[1])
rna_dna_m_c = UtC - rna_m_c - dna_m_c
#print(rna_m_c, dna_m_c, rna_dna_m_c)#5 6 4 

with open(logfile, 'a') as log:
        log.write("\ndist_func_multi: uniform\n")
        log.write("\nNumber of RNA multi counts is {}".format(rna_m_c))
        log.write("\nNumber of DNA multi counts is {}".format(dna_m_c))
        log.write("\nNumber of DNA and RNA multi counts is {}".format(rna_dna_m_c))
        
multi_countDNA = [ConstMulti_for_DNAmulti for i in range(dna_m_c)]
multi_countRNA = [ConstMulti_for_RNAmulti for i in range(rna_m_c)]
multi_countRNADNA = [ConstMulti_for_DNAmultiRNAmulti for i in range(rna_dna_m_c)]

uni_trueDNA = uni_true[:dna_m_c]
uni_trueRNA = uni_true[dna_m_c:(rna_m_c+dna_m_c)]
uni_trueRNADNA = uni_true[(rna_m_c+dna_m_c):]

with open(output1, 'w') as mtff:
    multi_list_dist = []
    multi_list_dist_wtrue = []
    #make multies for DNA part
    for true, count in zip(uni_trueDNA, multi_countDNA):
        ids_true = str(true[0])+"_"+str(true[1])+"_dna"
        mtff.write("\t".join([ids_true, str(true[0]), str(true[1])])+"\n")
        multi_list_dist_wtrue.append(abs(true[0]-true[1]))
        dna_coord = random.uniform(low=1, high=CL, size=count).round(decimals=0)
        dna_coord1 = [int(i) for i in dna_coord]
        for dna_false in dna_coord1:
            mtff.write("\t".join([ids_true, str(true[0]), str(dna_false)])+"\n")
            multi_list_dist_wtrue.append(abs(true[0]-dna_false))
            multi_list_dist.append(abs(true[0]-dna_false))
    #make multies for RNA part
    for true, count in zip(uni_trueRNA, multi_countRNA):
        ids_true = str(true[0])+"_"+str(true[1])+"_rna"#ids rna_dna
        mtff.write("\t".join([ids_true, str(true[0]), str(true[1])])+"\n")
        multi_list_dist_wtrue.append(abs(true[0]-true[1]))
        rna_coord = random.uniform(low=1, high=CL, size=count).round(decimals=0)
        rna_coord1 = [int(i) for i in rna_coord]
        for rna_false in rna_coord1:
            mtff.write("\t".join([ids_true, str(rna_false), str(true[1])])+"\n")
            multi_list_dist_wtrue.append(abs(true[1]-rna_false))
            multi_list_dist.append(abs(true[1]-rna_false))
    #make double multies for multiRNA part and multiDNApart
    for true, count in zip(uni_trueRNADNA, multi_countRNADNA):
        ids_true = str(true[0])+"_"+str(true[1])+"_rnadna"#ids rna_dna
        multi_list_dist_wtrue.append(abs(true[0]-true[1]))
        rna_coord = random.uniform(low=1, high=CL, size=count).round(decimals=0)
        rna_coord1 = [int(i) for i in rna_coord]
        rna_coord1.append(true[0])
        dna_coord = random.uniform(low=1, high=CL, size=count).round(decimals=0)
        dna_coord1 = [int(i) for i in dna_coord]
        dna_coord1.append(true[1])
        iter_no = count*count
        it = 0
        for rna_false in rna_coord1:#может быть и не фолз
            for dna_false in dna_coord1:
                it+=1
                mtff.write("\t".join([ids_true, str(rna_false), str(dna_false)])+"\n")
                multi_list_dist_wtrue.append(abs(dna_false-rna_false))
                if it != iter_no-1:
                    multi_list_dist.append(abs(dna_false-rna_false))
    

No_Mtf = len(multi_list_dist_wtrue)
with open(logfile, 'a') as log:
        log.write("\n\nNumber of multi coord cont {}\n".format(No_Mtf))

random_prob = ((rna_m_c+dna_m_c+rna_dna_m_c)/(rna_m_c*(ConstMulti_for_RNAmulti+1) + (dna_m_c*(ConstMulti_for_DNAmulti+1))+(rna_dna_m_c*(ConstMulti_for_DNAmultiRNAmulti+1)**2)))*100
with open(logfile, 'a') as log:
        log.write("\n\nProbability of random choice (before filter) {}\n".format(random_prob))
        

###visualisation
global BL
BL = 500000
global BC
BC = CL//BL + 1

plt.hist(multi_list_dist, bins = BC)
plt.title('Multi_false{}'.format(MaMC))
plt.savefig(pic1)

plt.hist(multi_list_dist_wtrue, bins = BC)
plt.title('Multi_false_true{}'.format(MaMC))
plt.savefig(pic2)

#plt.hist([i+1 for i in multi_count])
#plt.title('Multi_counts{}'.format(MaMC))
#plt.savefig(pic3)