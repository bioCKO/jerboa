import os
import pandas as pd
import numpy as np
from Bio import SeqIO
import subprocess
import shlex
import sys


def anno_fake_provean(fl, save):
    dfs = []
    n = 0
    for f in fl:
        n += 1
        if n % 1000 == 0:
            print(n)
        df = pd.read_table(f, header=0)

        # read in the real provean var score and annotated it in to PROVEAN column
        df['PROVEAN'] = np.random.randn(df.shape[0]) * 5
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
                    species_var_list.append(','.join([pos, var['Ref'], var['Alt']]) + '\n')
                except TypeError:
                    print(pos, var['Ref'], var['Alt'])
                    print(gene_id)
                    raise
        except AttributeError:
            pos = str(species_var_df['Start'])
            if pos == '0':
                pos = '1'
            species_var_list.append(','.join([pos, species_var_df['Ref'], species_var_df['Alt']]) + '\n')

        sf = open(os.path.join(dir_path, 'provean', gene_id, '%s_%s.fa' % (gene_id, s)), 'w')
        sf.write('> %s\n' % species_protein_id)
        sf.writelines(species_protein)
        sf.close()
        vf = open(os.path.join(dir_path, 'provean', gene_id, '%s_%s.var' % (gene_id, s)), 'w')
        vf.writelines(species_var_list)
        vf.close()
    return


def run_provean(gene, file_dir, submit, num_threads='12', overwrite=True):
    gene_dir = os.path.join(file_dir, gene)
    fl = [i for i in os.listdir(gene_dir) if gene in i]
    species = list(set([os.path.splitext(i)[0].split('_')[1] for i in fl if 'blast' not in i]))
    first = True
    log_f = open(os.path.join(gene_dir, 'provean_%s.log' % submit), 'w')
    for s in species:
        var_file = os.path.join(gene_dir, [i for i in fl if (s in i) and ('.var' in i)][0])
        fa_file = os.path.join(gene_dir, [i for i in fl if (s in i) and ('.fa' in i)][0])
        if os.path.exists(os.path.splitext(fa_file)[0] + '.provean') and not overwrite:
            continue
        save_fasta_set = os.path.join(gene_dir, gene + '.blast_set')
        if first:
            first = False
            # do first provean with blast, save fasta
            try:
                cmd = ['/home/hanqing/provean-1.1.5/src/provean',
                       '-q', fa_file,
                       '-v', var_file,
                       '-d', '/home/hanqing/db/nr',
                       '--psiblast', '/home/hanqing/ncbi-blast-2.7.1+/bin/psiblast',
                       '--cdhit', '/home/hanqing/cdhit/cd-hit',
                       '--blastdbcmd', '/home/hanqing/ncbi-blast-2.7.1+/bin/blastdbcmd',
                       '--save_supporting_set', save_fasta_set,
                       '--num_threads', num_threads,
                       '--tmp_dir', '/home/hanqing/tmp/',
                       '--quiet']
                output = subprocess.check_output(cmd, timeout=300)
                subprocess.check_call(shlex.split("date '+%A %W %Y %X'"))
                subprocess.check_call(['echo', "Command Finish\t" + ' '.join(cmd)])
                f = open(os.path.splitext(fa_file)[0] + '.provean', 'wb')
                f.write(output)
                f.close()
            except TimeoutError:
                log_f.write('# Time out\n')
                log_f.write(' '.join(cmd)+'\n')
                continue
        else:
            # do other provean without blast
            try:
                cmd = ['/home/hanqing/provean-1.1.5/src/provean',
                       '-q', fa_file,
                       '-v', var_file,
                       '-d', '/home/hanqing/db/nr',
                       '--psiblast', '/home/hanqing/ncbi-blast-2.7.1+/bin/psiblast',
                       '--cdhit', '/home/hanqing/cdhit/cd-hit',
                       '--blastdbcmd', '/home/hanqing/ncbi-blast-2.7.1+/bin/blastdbcmd',
                       '--subject_sequences', save_fasta_set + '.fasta',
                       '--tmp_dir', '/home/hanqing/tmp/',
                       '--quiet']
                output = subprocess.check_output(cmd, timeout=300)
                subprocess.check_call(shlex.split("date '+%A %W %Y %X'"))
                subprocess.check_call(['echo', "Command Finish\t" + ' '.join(cmd)])
                f = open(os.path.splitext(fa_file)[0] + '.provean', 'wb')
                f.write(output)
                f.close()
            except TimeoutError:
                log_f.write('# Time out\n')
                log_f.write(' '.join(cmd) + '\n')

    pfl = [i for i in os.listdir(gene_dir) if (gene in i) and ('provean' in i)]
    log_f.write('## Number of species: %d\n' % len(species))
    log_f.write('## Number of provean results: %d\n' % len(pfl))
    log_f.close()
    if len(species) != len(pfl):
        subprocess.check_call(['echo', "!Error\t" + "Gene: %s\n" % gene])
    else:
        subprocess.check_call(["echo", "!Success\t" + "Gene: %s\n" % gene])
    return


if __name__ == '__main__':
    submit = sys.argv[2]
    start = int(sys.argv[1]) * (int(sys.argv[2]) - 1)
    end = start + int(sys.argv[1])
    gene_list = open('/oasis/projects/nsf/csd473/hanqing/data/gene_list').readlines()
    gene_list = [i.strip() for i in gene_list]
    try:
        gene_list = gene_list[start:end]
    except IndexError:
        gene_list = gene_list[start:]

    for gene in gene_list:
        run_provean(gene=gene, file_dir='/oasis/projects/nsf/csd473/hanqing/data/jerboa/provean',
                    submit=submit, overwrite=False)

