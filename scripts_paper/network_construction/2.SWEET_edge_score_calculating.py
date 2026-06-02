import argparse
import numpy as np
import sys


def check_file(expres):
    checkset = set(["", "NA", "Na", "na", "nan", "null"])
    for c in checkset:
        loc = np.where(expres == c)
        if loc[0].size:
            expres[loc] = "0"
            print(f"Notice! There is {c} in the 'gene expression matrix' file and it will be assigned to 0.")
    return expres


parser = argparse.ArgumentParser(description="Manual")
parser.add_argument("-f", type=str, default="./example/expression.txt",
                    help="A path to 'gene expression matrix' file")  # gene expression matrix
parser.add_argument("-w", type=str, default="./example/weight.txt",
                    help="A path to 'sample weight' file (i.e., the output file from step 1)")  # sample weight file
parser.add_argument("-p", type=str, default="./example/patient.txt",
                    help="A path to 'samples of interest' file")  # samples of interest
parser.add_argument("-g", type=str, default="./example/gene.txt",
                    help="A path to 'genes of interest' file")  # genes of interest
parser.add_argument("-s", type=str, default="./example",
                    help="A path to the output 'confidence scores of edges' files for each sample of interest")  # output path

args = parser.parse_args()
file_e, file_w = args.f, args.w
file_p, file_g = args.p, args.g
save_path = (args.s).rstrip('/')

# open 'samples of interest' file
patset = set()
with open(file_p, mode='r') as rline:
    for nline in rline:
        tem = nline.strip('\n').split('\t')
        patset.add(tem[0])
del rline, nline, tem
if not patset:
    print("Warning! There is no sample ID in the 'samples of interest' file.")
    sys.exit()

# open 'genes of interest' file
geneset = set()
with open(file_g, mode='r') as rline:
    for nline in rline:
        tem = nline.strip('\n').split('\t')
        geneset.add(tem[0])
del rline, nline, tem
if not geneset:
    print("Warning! There is no gene ID in the 'genes of interest' file.")
    sys.exit()

# open 'sample weight' file
weight = {}
with open(file_w, mode='r') as rline:
    _ = rline.readline()
    for nline in rline:
        p, w, *_ = nline.strip('\n').split('\t')
        weight.update({p: float(w)})
del rline, nline, p, w, _
if not weight:
    print("Warning! There is no sample ID in the 'sample weight' file.")
    sys.exit()

# open 'gene expression matrix' file
gene, value = [], []
with open(file_e, mode='r') as rline:
    pat = rline.readline().strip('\n').split('\t')[1:]
    for nline in rline:
        g, *v = nline.strip('\n').split('\t')
        if g in geneset:
            value += v
            gene.append(g)
del rline, nline, g, v
patlen, genelen = len(pat), len(gene)
if (not patlen):
    print("Warning! The 'gene expression matrix' file is empty")
    sys.exit()

# check the 'samples of interest' and 'genes of interest' in expression file
patloc = [l for l, p in enumerate(pat) if p in patset]
if (not genelen) or (len(patloc) != len(patset)):
    print("Warning! The expression file cannot be mapped to 'samples of interest' or 'genes of interest' file")
    sys.exit()
if len(set(pat) & weight.keys()) != patlen:
    print("Warning! The sample ID(s) in the expression file cannot be mapped to 'sample weight' file")
    sys.exit()
del patset, geneset
print(f"patient : {len(patloc)}\ngene : {genelen}")

gene = np.array(gene)
value = np.array(value).reshape(genelen, patlen)
value = check_file(value)
value = value.astype(float)
loc = np.where(np.sum(value, axis=1) == 0)
if len(loc[0]) != 0:
    tem = ','.join(str(i)for i in gene[loc])
    print('Processing: delete gene(s) with zero expression values in all samples:'+tem)
    value = np.delete(value, loc, 0)
    gene = np.delete(gene, loc)

# value = np.array(value,dtype=float)


agg = np.corrcoef(value)

for l in patloc:
    p = pat[l]
    value_s = np.c_[value, value[:, l]]
    value_s = np.corrcoef(value_s)
    value_s = weight[p] * (value_s - agg) + agg
    with open(f"{save_path}/{p}.txt", mode='w') as wline:
        wline.write("gene1\tgene2\traw_edge_score\n")
        for l, g1, v1 in zip(range(genelen), gene, value_s):
            wline.write('\n'.join(g1+'\t'+g2+'\t'+str(v2)
                        for g2, v2 in zip(gene[(l+1):], v1[(l+1):])))
            wline.write('\n')

del value, agg, patloc, value_s, p, wline, weight, l, g1, v1

print("Finish")
