import pandas as pd
import csv
import numpy as np
from collections import Counter


# Read in mRNA file whole_mRNA_944samples.csv
def mRNA_Read(filepath):
    fea = []
    with open(filepath) as csvfile:
        csv_reader = csv.reader(csvfile)
        data_header = next(csv_reader)[0]
        sub = [x[5:] for x in data_header.split('\t')[1:]]
        for row in csv_reader:
            tmp = row[0].split('\t')[1:]
            fea.append(tmp)

    fea = np.transpose(np.array(fea, dtype=np.float32))
    return fea, sub


# Read in WGBS file whole_WGBS_Gene_985samples.csv
def wgbs_Read(filepath):
    fea = []
    with open(filepath) as csvfile:
        csv_reader = csv.reader(csvfile)
        data_header = next(csv_reader)[0]
        sub = [x for x in data_header.split(' ')[3:]]
        for row in csv_reader:
            tmp = row[0].split(' ')[3:]
            fea.append(tmp)

    fea = np.transpose(np.array(fea, dtype=np.float32))
    return fea, sub


# Read in miRNA file whole_miRNA_956sample.csv
def miRNA_Read(filepath):
    data = pd.read_csv(filepath, sep=',').values

    sub = []
    fea = []
    for i in range(len(data)):
        tmp1 = data[i][0].split()
        sub.append(tmp1[0])
        fea.append([float(x) for x in tmp1[1:]])
    fea = np.array(fea)
    return fea, sub


# Read excel, take in Subject_ID and Hip_Zscore to df
df = pd.read_excel('MM_Hip_Zscore2.xlsx')
sub_excel = df['Subject_ID'].tolist()
sub_excel_set = set(df['Subject_ID'].tolist())
lbl = df['Hip_Zscore'].to_numpy()

# Run defined functions
fea_mRNA, sub_mRNA = mRNA_Read('BoneData/whole_mRNA_944samples.csv')
fea_WGBS, sub_WGBS = wgbs_Read('BoneData/whole_WGBS_Gene_985samples.csv')
fea_mi, sub_mi = miRNA_Read('BoneData/whole_miRNA_956sample.csv')

# Create sets for subjects
sub_csv1 = set(sub_mRNA)
sub_csv2 = set(sub_WGBS)
sub_csv3 = set(sub_mi)

# Find common subjects in at least two files
common_subs = list(set().union(sub_csv1, sub_csv2, sub_csv3, sub_excel_set))

List_commonsubs = sub_mRNA + sub_WGBS + sub_mi + sub_excel
totalsub = len(sub_mRNA+sub_WGBS+sub_mi+sub_excel)

countCommon= dict(Counter(List_commonsubs))

duplicates = {key:value for key, value in countCommon.items() if value > 1}
print(duplicates)

dupCounter = 0
for key, value in duplicates.items():
    dupCounter = dupCounter + 1

with open('commonsubjects.txt', 'w') as f:
    f.write("Total Subjects: ")
    f.write(str(totalsub))
    f.write("\nCommon Subjects: ")
    f.write(str(dupCounter))
    f.write("\nCommon subjects and occurrences:\n")
    for key,value in duplicates.items():
        f.write(f"{key} {value}\n")

# read z score
df = pd.read_excel('MM_Hip_Zscore2.xlsx')
sub = df['Subject_ID'].tolist()
lbl = df['Hip_Zscore'].to_numpy()

highScore = []
lowScore = []
countHigh = 0
countLow = 0

# Find high and low scores
for i in range(len(sub)):
    if lbl[i] >= 0.8:
        highScore.append((sub[i], lbl[i]))
        countHigh = countHigh + 1
    else:
        lowScore.append((sub[i], lbl[i]))
        countLow = countLow + 1

tot = countHigh + countLow

# Sort scores by descending hip score
highScore = sorted(highScore, key=lambda x: x[1], reverse=True)
lowScore = sorted(lowScore, key=lambda x: x[1], reverse=True)

# Write results to highlowBMDScore.txt
with open('highlowBMDscores.txt', 'w') as f:

    f.write("Total Hip Bone Mineral Density Samples: ")
    f.write(str(tot)+'\n')

    f.write("\nHigh Hip BMD Score Count: ")
    f.write(str(countHigh))

    f.write('\nHigh BMD Score:\n')
    for sub, zscore in highScore:
        f.write(f"{sub} {zscore}\n")

    f.write("\nLow Hip BMD Score Count: ")
    f.write(str(countLow))
    f.write('\nLow BMD Score:\n')
    for sub, zscore in lowScore:
        f.write(f"{sub} {zscore}\n")

# Print
print("Number of High Hip BMD Scores:", countHigh, "\nNumber of Low Hip BMD Scores:", countLow)

'''
# Print and save common subjects to commonsubjects.txt
with open('commonsubjects.txt', 'w') as f:
    f.write("Common subjects:\n")
    for sub in common_subs:
        f.write("\n" + sub + "\n")
        # Print subject and features
        f.write("Features:\n")
        files = []

        # Creating functions for each case
        def mRNA_sub():
            files.append("mRNA")
            idx_mRNA = sub_mRNA.index(sub)
            f.write("mRNA:\n")
            f.write(str(fea_mRNA[idx_mRNA]) + "\n")
            return idx_mRNA

        def WGBS_sub():
            files.append("WGBS")
            idx_WGBS = sub_WGBS.index(sub)
            f.write("WGBS:\n")
            f.write(str(fea_WGBS[idx_WGBS]) + "\n")
            return idx_WGBS

        def mi_sub():
            files.append("miRNA")
            idx_mi = sub_mi.index(sub)
            f.write("miRNA:\n")
            f.write(str(fea_mi[idx_mi]) + "\n")
            return idx_mi

        def excel_sub():
            files.append("Excel")
            idx_excel = sub_excel.index(sub)
            f.write("Z_Score:\n")
            f.write(str(lbl[idx_excel]) + "\n")
            return idx_excel

        # Identifying which files share the subject, case by case
        # Using if rather than elif to make sure each condition is evaluated

        # Function to identify if more than 2 files share subject

        if sub in sub_mRNA and sub_WGBS:
            output = mRNA_sub()
            output2 = WGBS_sub()
        else:
            print("Common subject not in mRNA and WGBS")

        if sub in sub_mRNA and sub_mi:
            output = mRNA_sub()
            output2 = mi_sub()
        else:
            print("Common subject not in mRNA and miRNA")

        if sub in sub_mRNA and sub_excel:
            output = mRNA_sub()
            output2 = excel_sub()
        else:
            print("Common subject not in mRNA and Excel")

        if sub in sub_excel and sub_WGBS:
            output = excel_sub()
            output2 = WGBS_sub()
        else:
            print("Common subject not in Excel and WGBS")

        if sub in sub_excel and sub_mi:
            output = excel_sub()
            output2 = mi_sub()
        else:
            print("Common subject not in Excel and miRNA")

        if sub in sub_WGBS and sub_mi:
            output = WGBS_sub()
            output2 = mi_sub()
        else:
            print("Common subject not in WBGS and miRNA")

        if sub in sub_mRNA and sub_mi and sub_excel:
            output = mRNA_sub()
            output2 = mi_sub()
            output3 = excel_sub()
        else:
            print("Common subject not in mRNA, WGBS, Excel")

        if sub in sub_mRNA and sub_mi and sub_WGBS:
            output = mRNA_sub()
            output2 = mi_sub()
            output3 = WGBS_sub()
        else:
            print("Common subject not in mRNA, miRNA, WGBS")

        if sub in sub_mRNA and sub_excel and sub_WGBS:
            output = mRNA_sub()
            output2 = excel_sub()
            output3 = WGBS_sub()
        else:
            print("Common subject not in mRNA, Excel, WGBS")

        if sub in sub_mi and sub_excel and sub_WGBS:
            output = mi_sub()
            output2 = excel_sub()
            output3 = WGBS_sub()
        else:
            print("Common subject not in miRNA, Excel, WGBS")

        if sub in sub_mRNA and sub_mi and sub_excel and sub_WGBS:
            output = mRNA_sub()
            output2 = mi_sub()
            output3 = excel_sub()
            output4 = WGBS_sub()
        else:
            print("Common subject not in mRNA, miRNA, Excel, WGBS")'''

''' if sub in sub_WGBS:
            files.append("WGBS")
            idx_WGBS = sub_WGBS.index(sub)
            f.write("WGBS Features:\n")
            f.write("Data: " + str(fea_WGBS[idx_WGBS]) + "\n")
        if sub in sub_mi:
            files.append("miRNA")
            idx_mi = sub_mi.index(sub)
            f.write("miRNA Features:\n")
            f.write("Data: " + str(fea_mi[idx_mi]) + "\n")
        if sub in sub_excel:
            files.append("Excel")
            idx_excel = sub_excel.index(sub)
            f.write("Excel Features:\n")
            f.write("Data: " + str(lbl[idx_excel]) + "\n")

        f.write("Shared files: " + ", ".join(files) + "\n")'''