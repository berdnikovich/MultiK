import sys

up_prior = sys.argv[1]
down_prior = sys.argv[2]
result = sys.argv[3]

with open(down_prior, 'r') as p1, open(up_prior, 'r') as p2:
    with open(result, 'w') as out:
        upstr_list = []#dist - отрицательные, файлы rna-dna
        sum_prob = 0
        line_1 = p2.readline()
        up_d = float(line_1.strip().split()[1])
        line1 = p1.readline()
        down_d = float(line1.strip().split()[1])
        main_d = (up_d+down_d)/2
        sum_prob += main_d
        for i in p2:#up rna>dna
            #print(i)
            i_ = i.strip().split()
            inew = ["-"+str(i_[0]), float(i_[1])]
            upstr_list.append(inew)
            sum_prob+=inew[1]#здесь добавлено 0 значение
            #out.write("-" + i)
        upstr_list.reverse()
        upstrdown_list = []
        for ll in upstr_list:
            upstrdown_list.append(ll)
        upstrdown_list.append([0, main_d])
        #out.write("\t".join(["0", "0.000230644505309"])+"\n")
        for j in p1:
            #print(j)
            j_=j.strip().split()
            jnew = [j_[0], float(j_[1])]
            upstrdown_list.append(jnew)
            sum_prob+=jnew[1]
        all_list = [[str(f[0]), str(f[1]/sum_prob)] for f in upstrdown_list]#здесь уже нормированная вероятность 
        for line in all_list:
            out.write("\t".join(line)+"\n")
        print("Overall sum of prob in prior distrib is "+str(sum_prob))