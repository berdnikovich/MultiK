import bisect
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


outputPath = sys.argv[1]#"/home/berdnikovich/try_mHiC/result_simulation_simmple/new_generation/sim1"#sys.argv[1]
BL = int(sys.argv[2])
CL = int(sys.argv[3])

global BC
BC = int((CL+ BL) /  BL)
prior_up = os.path.join(outputPath, 'up_prior/s5_prior.mhic')
prior_down = os.path.join(outputPath, 'down_prior/s5_prior.mhic')
all_prior =  os.path.join(outputPath, 's5_prior12')
multi_file = os.path.join(outputPath, 's6/outFile.multi.mHiC')
uni_reads = os.path.join(outputPath, 's4/UNI_contacts_ALL_{}.Red_c'.format(BL))

##visualisation/tables
global map_prior_pic
map_prior_pic = os.path.join(outputPath, "s6/prior.png")
global heatmap_multi_prob
heatmap_multi_prob = os.path.join(outputPath, 's6/heatmap_multi_prob.png')
global heatmap_multi_post_max
heatmap_multi_post_max = os.path.join(outputPath, "s6/heatmap_multi_post_max.png")
global heatmap_multi_post_max_uni
heatmap_multi_post_max_uni = os.path.join(outputPath, "s6/heatmap_multi_post_max_uni.png")
global table_results
table_results = os.path.join(outputPath, "s6/table_post_pr_true_pair.csv")
global logfile
logfile = os.path.join(outputPath, "s6/parse_results.log")

def read_spline_prior(splineFilePath):
    """
    read in genomic distance prior
    """
    spline = {}
    splineLine = []# задаем переменные для цикла далее

    splineFile = open(splineFilePath, "r")
    for line in splineFile:
        splineLine = line.rstrip().split()
        spline[int(splineLine[0])] = float(splineLine[1])#здесь можем взять отрицательные дистанции к примеру
    interProb = min(spline.values())/2#вот здесь нам бы посмотреть, что дальше эта переменная творит - кажется это вероятность двугой хромосомы
    splineFile.close()
    return (spline, interProb)
dict_prior, interpr = read_spline_prior(all_prior)

def get_spline_probability(spline, interProb, pairString):    
    """
    init pi based on chr and pos distance
    """
    pair = [] 
    splineList = []
    
    pair = pairString.rstrip().split("_")
    ch1 = pair[0]
    pos1 = int(float(pair[1]))
    ch2 = pair[2]
    pos2 = int(float(pair[3]))
    #dist = abs(pos2 - pos1)
    dist = pos2 - pos1#pos1 = rna, pos2 = dna, если rna > dna - upstr, dist  будет отрицательным
    
    
    if ch1 != ch2:
        return interProb#кажется это деленная на 2 наименьшая вероятность
    elif dist in spline.keys():
        return spline[dist]
    else:#если дистанции нет в spline словаре
        splineList = list(spline.keys())
        indz = splineList.index(0)#вызываем индекс 0
        if dist < 0:#upstr
            #ind of dists is [0,indz] slice is [:indz] должен быть отсортирован, чтобы находился правильный индекс 
            splineLen1 = len(splineList[:indz])-1
            xPoint1 = min(bisect.bisect_right(splineList[:indz], dist), splineLen1 - 1)#ищем индекс сревниваемого 
            #элемента, если расстояние слишком большое берем наибольшее расстояние из списка 
            return spline[splineList[xPoint1]]
        else:#downstr
            #ind of dists is [indz+1] slice is [indz+1:] 
            splineLen2 = len(splineList[indz+1:])-1
            xPoint2 = min(bisect.bisect_right(splineList[indz+1:], dist), 0)#splineLen2-1
            return spline[splineList[xPoint2]]

def build_map_prior(CL1, BL1):#grid from CL, BL, 
    #counts sum in cell, BL in counts df should be less than BL
    #["st"]-rna part, ["end"]-dna part, "count"-count in counts_df
    no_of_bins = int((CL1+ BL1) /  BL1)
    grid_map = [i for i in range(1, no_of_bins+1)]
    map_df = pd.DataFrame(0, index = grid_map, columns = grid_map)
    for rna_ind in grid_map:
        for dna_ind in grid_map:
            a = "_".join(["chr20", str(rna_ind*BL1-BL1//2),"chr20", str(dna_ind*BL1-BL1//2)])
            map_df.loc[rna_ind, dna_ind] = get_spline_probability(dict_prior, interpr, a)#col is rna part, row is dna part
    return map_df
map_prior = build_map_prior(CL, BL)

fig1, mprior = plt.subplots(nrows=1, ncols=1)
mprior = sns.heatmap(map_prior)
mprior.set(xlabel="DNA-part (bin index)", ylabel="RNA-part (bin index)")
mprior.set_title("Prior probability", fontsize =20)
fig1.savefig(map_prior_pic, bbox_inches = 'tight')
plt.close(fig1)

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

df_multi = pd.read_csv(multi_file, sep='\t', names = ['id', 'chr1', 'st', 'chr2', 'end', 'count'])
#map_multi = build_map(CL, BL, df_multi)

#sns.heatmap(map_multi)
#plt.title("Probabilities for Multi reads", fontsize =20)
#plt.xlabel("DNA-part (bin index)")
#plt.ylabel("RNA-part (bin index)")
#plt.savefig(heatmap_multi_prob, bbox_inches = 'tight')

df_multi['max_count'] = df_multi.groupby(["id"])['count'].transform(max)
df_multi["max_prob"] = df_multi["max_count"] == df_multi["count"]
del df_multi['max_count']
df_multi["post_count"] = df_multi.apply(lambda row: 0 if row.max_prob == False else 1, axis=1)
del df_multi["max_prob"]

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
        count1 = row["post_count"]
        map_df.loc[rna_p, dna_p]= map_df.loc[rna_p, dna_p] + count1#col is rna part, row is dna part
    return map_df

map_multi_post = build_map(CL, BL, df_multi)

fig2, multi_max = plt.subplots(nrows=1, ncols=1)
multi_max = sns.heatmap(map_multi_post, vmin=0, vmax=10)
multi_max.set(xlabel="DNA-part (bin index)", ylabel="RNA-part (bin index)")
multi_max.set_title("Multi reads with max probabilities", fontsize =20)
fig2.savefig(heatmap_multi_post_max, bbox_inches = 'tight')
plt.close(fig2)

#df_uni = pd.read_csv(uni_reads, sep='\t', names = ['chr1', 'st', 'chr2', 'end', 'count'])#uni_reads

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
    map_df = map_multi_post
    for index, row in counts_df.iterrows():
        rna_p = search_bin1(row.st) 
        dna_p = search_bin1(row.end)
        count1 = row["count"]
        map_df.loc[rna_p, dna_p]= map_df.loc[rna_p, dna_p] + count1#col is rna part, row is dna part
    return map_df

#map_multi_post_uni = build_map(CL, BL, df_uni)
#sns.heatmap(map_multi_post_uni)
#plt.title("Multi reads and Uni reads", fontsize =20)
#plt.xlabel("DNA-part (bin index)")
#plt.ylabel("RNA-part (bin index)")
#plt.savefig(heatmap_multi_post_max_uni, bbox_inches = 'tight')

df_multi["true_pair"] = df_multi.apply(lambda row: 1 if row.id.split('_')[:2] == [str(int((row.st+BL//2)/BL)), str(int((row.end+BL//2)/BL))] else 0, axis=1)

df_multi = df_multi.sort_values(['id', 'count'], ascending=False)
df_multi['rank_prob'] = df_multi.groupby(['id'])['count'].rank('min', ascending=False)

df_multi.to_csv(table_results ,header=False)

df_model_true = df_multi.loc[(df_multi["post_count"] == 1) & (df_multi["true_pair"] == 1)]
result_prob = (len(df_model_true)/len(df_multi))*100
with open(logfile, 'w') as log:
        log.write("\n Prob for model: {}\n".format(result_prob))
