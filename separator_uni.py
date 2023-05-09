#dna-rna:up, rna-dna:down
import sys
import os
global BL
BL = 500000
outputPath = sys.argv[1]

global input1
input1 = os.path.join(outputPath, "UNI_contacts_ALL_{}.Red_c".format(BL))
global down_out
down_out = os.path.join(outputPath, "Uni_down_{}.Red_c".format(BL))
global up_out
up_out = os.path.join(outputPath, "Uni_up_{}.Red_c".format(BL))


#fi_uni = 'UNI_contacts_ALL_500000.Red_c'
#fi_uni_up_dna_rna = 'Uni_up_500000.Red_c'
#fi_uni_down_rna_dna = 'Uni_down_500000.Red_c'

with open(input1, 'r') as uni:
    with open(up_out, 'w') as up, open(down_out, 'w') as down:
        for i in uni:
            i_=i.strip().split('\t')
            rna = int(i_[1])
            dna = int(i_[3])
            if rna >= dna:#
                up.write("\t".join([str(i_[0]), str(i_[3]), str(i_[2]), str(i_[1]), str(i_[4])])+"\n")
            if rna <= dna:
                down.write(i)
