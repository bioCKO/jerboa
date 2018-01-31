import os
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import MutableSeq


def anno_fake_provean(fl, save):
    dfs = []
    n = 0
    for f in fl:
        n += 1
        if n % 1000 == 0:
            print(n)
        df = pd.read_table(f, header=0)

        # read in the real provean var score and annotated it in to PROVEAN column
        df['PROVEAN'] = np.random.randn(df.shape[0])*5
        dfs.append(df)

    tdf = pd.concat(dfs)
    tdf.to_csv(save, index=False, sep='\t',
               columns=['Jaculus_jaculus protein', 'Start', 'End',
                        'Ref', 'End.1', 'Identity', 'Orthologous protein', 'PROVEAN'])


def get_provean_input(gene_id, dir_path):
    """
    get the modified fasta seq
    get the var list
    save as separate file

    :param gene_id:
    :return:
    """
    protein_fa_dict = SeqIO.to_dict(SeqIO.parse(
        os.path.join(dir_path, 'pep', '%s.fa' % gene_id), format='fasta'),
        key_function=lambda i: i.id.split('.')[0])

    var_list = pd.read_table(os.path.join(dir_path, 'var_ortho', '%s.var.tsv' % gene_id), header=0,
                             index_col='Orthologous species', keep_default_na=False, na_values='??')

    species = var_list.index.unique().tolist()
    os.mkdir(os.path.join(dir_path, 'provean', gene_id))
    for s in species:
        species_var_df = var_list.loc[s]
        try:
            species_protein_id = species_var_df['Orthologous protein'].iloc[0]
            species_protein = protein_fa_dict[species_protein_id]
        except AttributeError:
            if isinstance(species_var_df, pd.Series):
                species_protein_id = species_var_df['Orthologous protein']
                species_protein = protein_fa_dict[species_protein_id]
            else:
                raise
        species_var_list = []
        try:
            for idx, var in species_var_df.iterrows():
                pos = str(var['Start'])
                if pos == '0':
                    pos = '1'
                try:
                    species_var_list.append(','.join([pos, var['Ref'], var['Alt']])+'\n')
                except TypeError:
                    print(pos, var['Ref'], var['Alt'])
                    print(gene_id)
                    raise
        except AttributeError:
            pos = str(species_var_df['Start'])
            if pos == '0':
                pos = '1'
            species_var_list.append(','.join([pos, species_var_df['Ref'], species_var_df['Alt']])+'\n')

        sf = open(os.path.join(dir_path, 'provean', gene_id, '%s_%s.fa' % (gene_id, s)), 'w')
        sf.write('> %s\n' % species_protein_id)
        sf.writelines(species_protein)
        sf.close()
        vf = open(os.path.join(dir_path, 'provean', gene_id, '%s_%s.var' % (gene_id, s)), 'w')
        vf.writelines(species_var_list)
        vf.close()
    return


fl = [i for i in os.listdir('/Users/hq/JacJac/var_ortho') if 'var' in i]
n = 0
for f in fl:
    n += 1
    if n % 1000 == 0:
        print(n)

    gid = f.split('.')[0]
    get_provean_input(gid, dir_path='/Users/hq/JacJac/')

