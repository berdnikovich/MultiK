import sys
import os
import pandas as pd
from collections import Counter
import numpy
import seaborn as sns
import matplotlib.pyplot as plt

#global CL
#CL = 64444167
#global BL
#BL = 500000

outputPath = sys.argv[1]#'/home/berdnikovich/try_mHiC/simulation/10_sim/sim1/s4
picPath = sys.argv[2]#'/home/berdnikovich/try_mHiC/simulation/10_sim/sim1/visualisation
CL = int(sys.argv[3])
BL = int(sys.argv[4])
libName = '3stage_{}'.format(BL)
global logfile
logfile = os.path.join(outputPath, libName+".fithic.log")

global BC
BC = int((CL+ BL) /  BL)

###filenames
global frfile
frfile = os.path.join(outputPath, "fragments_{}".format(BL))
global unifile
unifile = os.path.join(outputPath, 'Uni_rna_dna_coo')
global output1
output1 = os.path.join(outputPath, "UNI_contacts_ALL_{}.Red_c".format(BL))
global multifile
multifile = os.path.join(outputPath, 'Multi_tf')
global output2
output2 = os.path.join(outputPath, "MULTI_contacts_ALL_{}.Red_c".format(BL))
##visualisation/tables
global unitab
unitab = os.path.join(picPath, "UNI_contacts_ALL_table_{}".format(BL))
global unipic
unipic = os.path.join(picPath, 'UNI_heatmap_{}.png'.format(BL))
global multitab
multitab = os.path.join(picPath, "MULTI_contacts_ALL_table_{}".format(BL))
global multipic
multipic = os.path.join(picPath, 'Multi_heatmap_{}.png'.format(BL))
global multipic_true
multipic_true = os.path.join(picPath, "true_multi_map.png")
global multipic_false
multipic_false = os.path.join(picPath, "false_multi_map.png")



with open(frfile, 'w') as fr_file:
    curr_chr = 'chr20'
    if ((CL % BL) == 0):
        interval_end = CL
    else:
        interval_end = (int((CL+ BL) /  BL)) *  BL
        for val in range(0, interval_end,  BL):
            curr_start = val
            curr_mid = val +  BL / 2
            if val > 0:
                fr_file.write('\n')
            fr_file.write(str(curr_chr) +'\t'+str(curr_start)+'\t'+str(int(curr_mid))+'\t'+str(1)+'\t'+str(1)) 
            
bin_line = [(i*BL, (i+1)*BL) for i in range(BC)]
def search_bin(coo):
    for bi_i, bi in enumerate(bin_line):
        if coo > bi[0] and coo <= bi[1]:#включаем правый край
            bin_indx = bi_i
            break
    return bin_indx+1#индексы в бинах начинаются с 1

##uni-reads

with open(unifile, 'r') as unif:#bin ind from 1
    NoUni = 0
    No_bb = 0#number of uni contacts w/in one bin
    li_uni_bins = []
    for uni_coo in unif:
        uni_coo_ = uni_coo.strip().split("\t")
        rna_coo = int(uni_coo_[0])
        dna_coo = int(uni_coo_[1])
        rna_bin = search_bin(rna_coo)#bin ind from 1
        dna_bin = search_bin(dna_coo)#bin ind from 1
        if rna_bin == dna_bin:
            No_bb+=1
            li_uni_bins.append([rna_bin, dna_bin])
        else:
            NoUni+=1
            li_uni_bins.append([rna_bin, dna_bin])

rna_bins_u = [i[0] for i in li_uni_bins]
dna_bins_u = [i[1] for i in li_uni_bins]
df_rna_dna = pd.DataFrame({'st' : rna_bins_u,'end' : dna_bins_u})
df_rna_dna['Count'] = 1
df_rna_dna1 = df_rna_dna.groupby(['st','end']).Count.count().reset_index()

li_ch = ['chr20' for i in range(len(df_rna_dna1))]
df_uni_rna_dna = pd.DataFrame({'ch1': li_ch, 'st' : df_rna_dna1['st']*BL-BL/2,# rna='start'
                                        'ch2': li_ch ,'end' : df_rna_dna1['end']*BL-BL/2, #dna='end'
                                        'count' : df_rna_dna1['Count']})

df_uni_rna_dna = df_uni_rna_dna.astype({'st': 'int32', 'end': 'int32'})
#for uni list use mid of bin
df_uni_rna_dna.to_csv(output1, sep='\t', header = False, index = False)

def build_map(CL1, BL1, counts_df):#grid from CL, BL, 
    #counts sum in cell, BL in counts df should be less than BL
    #["st"]-rna part, ["end"]-dna part, "count"-count in counts_df
    no_of_bins = int((CL1+ BL1) /  BL1)
    bin_line1 = [(i*BL1, (i+1)*BL1) for i in range(no_of_bins)]
    def search_bin1(coo):
        for bi_i, bi in enumerate(bin_line1):
            if coo > bi[0] and coo <= bi[1]:#включаем правый край
                bin_indx = bi_i
                break
        return bin_indx+1
    grid_map = [i for i in range(1, no_of_bins+1)]
    map_df = pd.DataFrame(0, index = grid_map, columns = grid_map)
    for index, row in counts_df.iterrows():
        rna_p = search_bin1(row.st) 
        dna_p = search_bin1(row.end)
        count1 = row["count"]
        map_df.loc[rna_p, dna_p]= map_df.loc[rna_p, dna_p] + count1#col is rna part, row is dna part
    return map_df

df_uni  = build_map(CL, BL, df_uni_rna_dna)

fig1, u = plt.subplots(nrows=1, ncols=1)
u = sns.heatmap(df_uni, vmin=0, vmax=10)
u.set(xlabel="DNA-part (bin index)", ylabel="RNA-part (bin index)")
u.set_title("UNIreads bin is {}kbp".format(BL//1000), fontsize =20)
fig1.savefig(unipic, bbox_inches = 'tight')
plt.close(fig1)


with open(logfile, 'w') as log:
        log.write("\n\nStart bin data\n")
        log.write("\n frag file: {}\n".format(logfile))
        log.write("\n rna_dna_cont uni: {}\n".format(NoUni))
        log.write("\n bin(rna) == bin(dna) uni: {}\n".format(No_bb))
##multi-reads

with open(multifile, 'r') as multif:#bin ind from 1
    No_deleted_mistake_in_true = 0
    li_multi_bins = []
    li_multi_true_bins = []
    li_multi_false_bins = []
    NoMulti = 0
    cank_no = 0
    count_true = 0
    count_true_false = 0
    pr_name = "_"
    chank_li = []
    #chank_li_li = []
    for multi_coo in multif:
        NoMulti+=1
        multi_coo_ = multi_coo.strip().split("\t")
        true_pair_new= multi_coo_[0].split("_")
        true_pair = "_".join(true_pair_new[:2])
        mmpart = true_pair_new[-1]
        if true_pair != pr_name:
            cank_no+=1
            pr_name = true_pair
            #chank_li_li.append(chank_li)
            chank_li = []
        true_pair_ = true_pair.split("_")
        true_rna = int(true_pair_[0])
        true_dna = int(true_pair_[1])
        new_name = str(search_bin(true_rna))+"_"+str(search_bin(true_dna))#bins from 1
        new_name1 = str(search_bin(true_rna))+"_"+str(search_bin(true_dna))+"_"+str(cank_no)+"_"+mmpart
        #пронумерованы чанки происходящие из одной рнк затравки
        name_pair = multi_coo_[1]+"_"+multi_coo_[2]#старые координаты
        rna_coo = int(multi_coo_[1])
        dna_coo = int(multi_coo_[2])
        rna_bin = search_bin(rna_coo)
        dna_bin = search_bin(dna_coo)
        new_name_pair = str(rna_bin)+"_"+str(dna_bin)
        if name_pair == true_pair:#тру пара в любом случае заносится
            li_multi_bins.append([new_name1,rna_bin, dna_bin])
            li_multi_true_bins.append([new_name1,rna_bin, dna_bin])
            count_true+=1
            count_true_false+=1
        else:#тру пара отработала и уже ушла - сюда стекаются все фолз пары
            if new_name_pair in chank_li:#бины у тру пары совпали с бинами у фолз пары
                No_deleted_mistake_in_true+=1#все фолз пары которые повторились с тру парой или любой другой парой
                continue# мы их отфильтровали, чтобы была тру пара только одна
            else:
                count_true_false+=1
                li_multi_bins.append([new_name1,rna_bin, dna_bin])
                li_multi_false_bins.append([new_name1,rna_bin, dna_bin])
        chank_li.append(new_name_pair)#попадает в тч тру пара
        
#отфильтровываются пары, которые повторяются между собой
#отфильтровываются пары, которые повторяют тру пару
#индексируются чанки

with open(logfile, 'a') as log:
        log.write("\n\nNumber of false pairs deleted due to duplication = {}\n".format(No_deleted_mistake_in_true))
        log.write("\n\nNumber of Multi pairs = {}\n".format(len(li_multi_bins)))
        log.write("\n\nRandom probabitity = {}\n".format(count_true/count_true_false*100))
li_id_multi = [i[0] for i in li_multi_bins]
li_rna_multi = [i[1] for i in li_multi_bins]
li_dna_multi = [i[2] for i in li_multi_bins]

li_ch = ['chr20' for i in range(len(li_dna_multi))]
df_rna_dna_m = pd.DataFrame({'id': li_id_multi, 'ch1': li_ch, 'st' : [int(i*BL-BL/2) for i in li_rna_multi], 'ch2': li_ch,'end' : [int(i*BL-BL/2) for i in li_dna_multi]})

df_rna_dna_m.to_csv(output2, sep='\t', header = False, index = False)

##visualisation/tables
df_rna_dna_m['count'] = 1
map_multi_df = build_map(CL, BL, df_rna_dna_m)

fig2, multi_alls = plt.subplots(nrows=1, ncols=1)
multi_alls = sns.heatmap(map_multi_df, vmin=0, vmax=5)
multi_alls.set(xlabel="DNA-part (bin index)", ylabel="RNA-part (bin index)")
multi_alls.set_title("Multireads bin is {}kbp".format(BL//1000), fontsize =20)
fig2.savefig(multipic, bbox_inches = 'tight')
plt.close(fig2)

li_id_multi_true = [i[0] for i in li_multi_true_bins]
li_rna_multi_true = [i[1] for i in li_multi_true_bins]
li_dna_multi_true = [i[2] for i in li_multi_true_bins]
li_ch = ['chr20' for i in range(len(li_multi_true_bins))]
df_rna_dna_m_true = pd.DataFrame({'id': li_id_multi_true, 'ch1': li_ch, 'st' : [int(i*BL-BL/2) for i in li_rna_multi_true], 'ch2': li_ch,'end' : [int(i*BL-BL/2) for i in li_dna_multi_true]})
df_rna_dna_m_true['count'] = 1

map_multi_true_df = build_map(CL, BL, df_rna_dna_m_true)

fig3, multi_trues = plt.subplots(nrows=1, ncols=1)
multi_trues = sns.heatmap(map_multi_true_df, vmin=0, vmax=5)
multi_trues.set(xlabel="DNA-part (bin index)", ylabel="RNA-part (bin index)")
multi_trues.set_title("True multireads bin is {}kbp".format(BL//1000), fontsize =20)
fig3.savefig(multipic_true, bbox_inches = 'tight')
plt.close(fig3)

li_id_multi_false = [i[0] for i in li_multi_false_bins]
li_rna_multi_false = [i[1] for i in li_multi_false_bins]
li_dna_multi_false = [i[2] for i in li_multi_false_bins]
li_ch = ['chr20' for i in range(len(li_multi_false_bins))]
df_rna_dna_m_false = pd.DataFrame({'id': li_id_multi_false, 'ch1': li_ch, 'st' : [int(i*BL-BL/2) for i in li_rna_multi_false], 'ch2': li_ch,'end' : [int(i*BL-BL/2)  for i in li_dna_multi_false]})
df_rna_dna_m_false['count'] = 1

map_multi_false_df = build_map(CL, BL, df_rna_dna_m_false)

fig4, multi_falses = plt.subplots(nrows=1, ncols=1)
multi_falses = sns.heatmap(map_multi_false_df, vmin=0, vmax=5)
multi_falses.set(xlabel="DNA-part (bin index)", ylabel="RNA-part (bin index)")
multi_falses.set_title("False multireads bin is {}kbp".format(BL//1000), fontsize =20)
fig4.savefig(multipic_false, bbox_inches = 'tight')
plt.close(fig4)