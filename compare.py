import os
from scipy import stats
import pandas as pd
from tqdm import tqdm
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import argparse
import io
def chrom_spliter(pred):
    res = {}
    res['I'] = pred[0:15072434]
    res['II'] = pred[15072434:30351855]
    res['III'] = pred[30351855:44135656]
    res['IV'] = pred[44135656:61629485]
    res['V'] = pred[61629485:82553665]
    res['X'] = pred[82553665:]
    return res


def compare(lst1, lst2, n,alpha):
    ref = sum(lst1) * n * alpha
    obs = sum(lst2)
    fc = obs/ref
    if obs <= n and ref/n < 1:
        p = stats.binomtest(obs, n, ref/n).pvalue
    else:
        p = 0
    return p, fc,obs,ref

def read_vcf(path_a):
    with open(path_a, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    )

def remove_background(file,bg_file):
    big = pd.concat([file,bg_file])
    big = big.groupby(['#CHROM', 'POS']).filter(lambda x: len(x) == 1)
    return big

def main(p_path, file_path, info_path, n, outpath,bg_path=False):
    # read model
    ref = chrom_spliter(pd.read_pickle(p_path)['mutation_factor'].to_list())
    obs = pd.DataFrame()
    file = os.listdir(file_path)
    info = pd.read_csv(info_path, names=['gene', 'chr', 'pos1', 'pos2'])
    if bg_path:
        bg = pd.read_csv(bg_path, skiprows=list(range(64)), sep='\t')
    for i in tqdm(file, desc="Removing Background"):
        if i.split('.')[-1] == 'vcf':
            try:
              a = read_vcf(file_path + "/{}".format(i))
              if bg_path:
                  try:
                      a = remove_background(a,bg)
                  except:
                      a = a
                      print(i)
              obs = pd.concat([obs, a])
            except:
              print(i)
    data = list(zip(obs['#CHROM'].to_list(), obs['POS'].to_list()))
    obs_dict = {}
    for k, v in ref.items():
        obs_dict[k] = [0] * len(v)
    # make result files
    for i in tqdm(data, desc="Creating mutation list"):
      try:
        obs_dict[i[0]][i[1] - 1] += 1
      except:
        continue

    for k in tqdm(obs_dict.keys(), desc="Preparing for annotation"):
        v = obs_dict[k]
        for i in range(len(v)):
            if v[i] > int(n):
                v[i] = 0
        obs_dict[k] = v
    sr = 0
    so = 0
    for k in tqdm(ref.keys(),desc="Adjusting mutation efficiency"):
        sr += sum(ref[k])
        so += sum(obs_dict[k])
    alpha = (so/len(file)) / sr
    print('alpha = {}'.format(alpha))
    pos = list(zip(info['pos1'].to_list(), info['pos2'].to_list()))
    chrom = info['chr'].to_list()
    gene = info['gene'].to_list()
    plst = []
    fclst = []
    genelst = []
    obslst = []
    reflst = []
    lgfclst = []
    for i in tqdm(range(len(pos)),desc="Comparing"):
        if chrom[i] in ['I','II','III','IV','V','X']:
          try:
            p, fc ,o,r= compare(ref[chrom[i]][pos[i][0]:pos[i][1]], obs_dict[chrom[i]][pos[i][0]:pos[i][1]], len(file),alpha)
            if o != 0:
                if fc != 0:
                    lgfc = math.log(fc,2)
                if p != 0:
                    p = -math.log(p,10)
                else:
                    p = 40
                plst.append(p)
                fclst.append(fc)
                obslst.append(o)
                reflst.append(r)
                lgfclst.append(lgfc)
                genelst.append(gene[i])
          except:
            print(gene[i])
    for i in range(len(plst)):
        if plst[i] > 40:
            plst[i] = 40
    plt.scatter(lgfclst,plst,s=5)
    plt.xlabel('Fold Change')
    plt.ylabel('P_val')
    plt.savefig(outpath+'/volcano.jpg')
    outable = pd.DataFrame({'gene':genelst,"observed":obslst,"expected":reflst,"foldchange":fclst,"p_value":plst})
    outable.to_csv(outpath+"/{}".format('volcano_table.csv'))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "EMS mutagenesis analyzation based on RF model")
    parser.add_argument("--model", "-m",help='model file')
    parser.add_argument("--data", "-d",help="vcf file folder")
    parser.add_argument("--ref","-r",help="gene range file")
    parser.add_argument("--threshold","-t",help="threshold when removing background")
    parser.add_argument("--out", "-o",help="output path")
    parser.add_argument("--background","-b",help="background path (optional)")
    args = parser.parse_args()
    main(args.model, args.data, args.ref,args.threshold,args.out,args.background)
