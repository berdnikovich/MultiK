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

global BL
BL = 500000
global CL
CL = 64444167
global BC
BC = CL//BL + 1
global NoR
NoR = 500
global PM
PM = 10
global ANoC
ANoC = 100
global exp_scale
lambd = -math.log(0.01)/(CL-BL)
exp_scale = 1/lambd 

outputPath = sys.argv[1]#'/home/berdnikovich/try_mHiC/simulation/10_sim/sim1/s4
libName = '1stage{}'.format(str(NoR))
picPath = sys.argv[2]#'/home/berdnikovich/try_mHiC/simulation/10_sim/sim1/visualisation

####
global logfile
logfile = os.path.join(outputPath, libName+".fithic.log")##вот так правильно в питоне пути создавать 
global shit_file1
shit_file1 = os.path.join(outputPath, "UNI_prime")
global out_file1
out_file1 = os.path.join(outputPath, "Uni_rna_dna_coo")
global out_file2
out_file2 = os.path.join(outputPath, "Multi_true_coo")
######
##visualisation
global pic1
pic1 = os.path.join(picPath, 'Stage1_uni_dist{}.png'.format(NoR))
global pic2
pic2 = os.path.join(picPath, 'Stage1_multi_true_dist{}.png'.format(NoR))


with open(logfile, 'w') as log:
        log.write("\n\nStart gene cont\n")
        log.write("\ndist_func: exp_scale = {}, loc = {}\n".format(exp_scale, '1'))
        log.write("\ncont_num_per_RNA: poisson, lam = {}\n".format(ANoC))
        
coo_rnas = random.uniform(low=1, high=CL, size=NoR).round(decimals=0)
counts_cont = random.poisson(lam=ANoC, size=NoR)

UNI_prime = open(shit_file1, 'w')
UNI_prime_li = []
NoC = 0
for bin1, M in zip(coo_rnas, counts_cont):
    i = 0# количество удачных совпадений
    while i < M:
        up_or_down = rn.choices(['Up', 'Down'], weights=[50, 50])# относительно bin1
        if up_or_down[0] == 'Up':
            len_of_contact = int(expon.rvs(scale=exp_scale,loc=1,size=1).round(decimals=0)[0])#расстояние может попасть в бин
            #print(len_of_contact)
            if bin1 - len_of_contact > 0:
                bin2 = bin1 - len_of_contact
                UNI_prime.write("\t".join([str(bin1), str(bin2)])+"\n")#dna-rna
                UNI_prime_li.append([str(int(bin1)), str(int(bin2))])
                NoC+=1
                i+=1
        else:
            len_of_contact = int(expon.rvs(scale=exp_scale,loc=1,size=1).round(decimals=0)[0])#расстояние может попасть в бин
            #print(len_of_contact)
            if bin1 + len_of_contact < CL:
                bin2 = bin1 + len_of_contact
                UNI_prime.write("\t".join([str(bin1), str(bin2)])+"\n")# rna-dna записываем в файл всегда от меньшего бина к большему
                UNI_prime_li.append([str(int(bin1)), str(int(bin2))])
                NoC+=1
                i+=1
UNI_prime.close()

random.shuffle(UNI_prime_li)
with open(logfile, 'a') as log:
        log.write("\n\nUni_prime = \n"+str(NoC))
NoMt =int(round(NoC*PM/100, 0))

Multi_true = UNI_prime_li[:NoMt]
Uni_li = UNI_prime_li[NoMt:]

with open(logfile, 'a') as log:
        log.write("\n\nMulti = \n"+str(len(Multi_true)))
        log.write("\n\nUni_end = \n"+str(len(Uni_li)))
        
with open(out_file1, 'w') as unif, open(out_file2, 'w') as multitf:
    for u in Uni_li:
        unif.write("\t".join(u)+"\n")
    for mt in Multi_true:
        multitf.write("\t".join(mt)+"\n")
        
##visualisation        


with open(out_file1, 'r') as unif, open(out_file2, 'r') as multitf:
    uni_dist_li = []
    multi_true_dist = []
    for u in unif:
        u_ = u.strip().split('\t')
        di = abs(int(u_[0])-int(u_[1]))
        uni_dist_li.append(di)
    for mt in multitf:
        mt_ = mt.strip().split('\t')
        di = abs(int(mt_[0])-int(mt_[1]))
        multi_true_dist.append(di)

plt.hist(uni_dist_li, bins = BC)
plt.title('uni_dist{}'.format(NoR))
plt.savefig(pic1)

plt.hist(multi_true_dist, bins = BC)
plt.title('multi_true_dist{}'.format(NoR))
plt.savefig(pic2)