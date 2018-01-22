import os
import pandas as pd
import numpy as np

dir_path = '/Volumes/Data/JacJac/var_fake'
fl = [os.path.join(dir_path, i) for i in os.listdir(dir_path) if 'var' in i]


def anno_fake_provean(fl):
    for f in fl:
        df = pd.read_table(f, header=0)

        # read in the real provean var score and annotated it in to PROVEAN column
        df['PROVEAN'] = np.random.randn(df.shape[0])*5

        df.to_csv(f[:-3]+'p.tsv', index=False, sep='\t')


fl = [os.path.join(dir_path, i) for i in os.listdir(dir_path) if 'var.p' in i]
def merge_var(fl):
    for f in fl:
        df = pd.read_csv(f, header=0, index_col='Orthologous species')
        species = df.index.unique()
        print(species)
        break
merge_var(fl)


