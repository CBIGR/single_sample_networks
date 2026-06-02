import argparse
import numpy as np
import sys
import os


def check_file(expres):
    checkset = set(["", "NA", "Na", "na", "nan", "null"])
    for c in checkset:
        loc = np.where(expres == c)
        if loc[0].size:
            expres[loc] = "0"
            print(f"There is {c} in the 'gene expression matrix' file and it will be assigned to 0.")
    return expres


parser = argparse.ArgumentParser(description="Manual")
parser.add_argument(
    "-p", type=str, default="./example/patient.txt", help="A path to 'samples of interest' file")
parser.add_argument("-l", type=str, default="./example", help="A path to the 'confidence scores of edges' file for each sample of interest (i.e., the output files from step 2)")
parser.add_argument(
    "-s", type=str, default="./example/mean_std.txt", help="A path to the output file(s)")
parser.add_argument("-z", type=bool, default=False,
                    help="Indicates whether the calculation of z score (Ture) or not (False)")

args = parser.parse_args()
file_p, file_l = args.p, (args.l).rstrip('/')
save = args.s
calculate_z = args.z

patlist = []
with open(file_p, mode='r') as rline:
    for nline in rline:
        tem = nline.strip('\n').split('\t')
        patlist.append(tem[0])

geneset = set()
pair = []
file = f"{file_l}/{patlist[0]}.txt"
if not os.path.exists(file):
    print(f"{file} not found")
    sys.exit()
with open(f"{file_l}/{patlist[0]}.txt", mode='r') as rline:
    _ = rline.readline()
    for nline in rline:
        if nline != '\n':
            val = nline.strip('\n').split('\t')
            geneset.add(val[0]+'\t'+val[1])
            pair.append(val[2])

for p in patlist[1:]:
    file = f"{file_l}/{p}.txt"
    if not os.path.exists(file):
        print(f"{file} not found")
        sys.exit()
    with open(f"{file_l}/{p}.txt", mode='r') as rline:
        _ = rline.readline()
        for nline in rline:
            if nline != '\n':
                val = nline.strip('\n').split('\t')
                if (val[0]+'\t'+val[1]) not in geneset:
                    print(f"Warning! In the sample {p}, there are gene pair(s) that cannot be found in the other samples.")
                pair.append(val[2])
pair = np.array(pair)
pair = check_file(pair)
pair = pair.astype(float)
vmean, vstd = np.mean(pair), np.std(pair)

with open(save, mode='w') as wline:
    wline.write(f"mean\t{vmean}\nstd\t{vstd}\n")

if calculate_z:
    for p in patlist:
        file = f"{file_l}/{p}.txt"
        if not os.path.exists(file):
            print(f"{file} not found")
            sys.exit()
        with open(f"{file_l}/{p}.txt", mode='r') as rline, open(f"{file_l}/{p}_zscore.txt", mode='w') as wline:
            _ = rline.readline()
            wline.write('gene1\tgene2\tz_score\n')
            for nline in rline:
                if nline != '\n':
                    val = nline.strip('\n').split('\t')
                    z = str(float(val[2])-vmean/vstd)
                    wline.write(f'{val[0]}\t{val[1]}\t{z}\n')

print("Finish")
