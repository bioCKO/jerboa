"""
1. Given a bunch of fasta files for several species, a table contains target genes and their orthologous gene IDs,
2. Find fasta sequence for each target gene's orthologous group,
   save sequences in a .fa file in the name of target gene ID,
3. Use MUSCLE to do MSA, save alignment result in a .aln file in the name of target gene ID,
4. Call variants from MSA follow the instruction of HGVS http://varnomen.hgvs.org/,
   Save a list of the variants in .var file in the name of target gene ID, add header indicates the species info,
   Save the variants in .vcf file
"""

import os
import subprocess
import multiprocessing
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio import AlignIO


def get_fasta_dict(fasta_dir):
    fasta_path_list = [os.path.join(fasta_dir, i) for i in os.listdir(fasta_dir) if 'fa' in i]
    species_list = [i.split('.')[0] for i in os.listdir(fasta_dir) if 'fa' in i]
    fasta_dict = {}
    for i in range(len(species_list)):
        fasta_dict[species_list[i]] = SeqIO.to_dict(SeqIO.parse(fasta_path_list[i], "fasta"),
                                                    key_function=lambda i: i.id.split('.')[0])
    return fasta_dict


def get_gene_ids(id_path, usecols, rename_col, query_col):
    """
    return a dict like this:
    {
        'GENE_ID':
            {
            'SPECIES': 'ortho gene ID' # including the target species and others
            ...
            }
        ...
    }

    :param id_path:
    :param usecols:
    :return:
    """
    df = pd.read_table(id_path, index_col='Gene stable ID', usecols=usecols)
    query_df = pd.read_table(id_path, index_col='Gene stable ID', usecols=query_col)
    query_ids = [list(set(i)) for i in query_df.values.tolist()]
    query_list = []
    for line in query_ids:
        if len(line) > 2:
            print(line)
        flag = True
        for j in line:
            if isinstance(j, str):
                query_list.append(j)
                flag = False
                break
        if flag:
            query_list.append(np.nan)

    df['Jaculus_jaculus protein'] = query_list

    col_dict = {}
    for i in range(len(usecols)):
        col_dict[usecols[i]] = rename_col[i]
    df.rename(columns=col_dict, inplace=True)

    _id_dict = df.to_dict('split')
    total_dict = {}
    for i in range(len(_id_dict['index'])):
        gene_id = _id_dict['index'][i]
        total_dict[gene_id] = {}
        for j in range(len(_id_dict['columns'])):
            total_dict[gene_id][_id_dict['columns'][j]] = _id_dict['data'][i][j]
    return total_dict


def get_orthogroup_seq(gene_id_dict, fasta_dict, save_path):
    species_list = list(fasta_dict.keys())
    not_in_idlist = [i for i in species_list if i + ' protein' not in list(gene_id_dict.values())[0]]
    for i in not_in_idlist:
        del species_list[species_list.index(i)]
    nf = open(os.path.join(save_path, 'error_id.txt'), 'w')
    counter = 0
    for gene_id, orthos in gene_id_dict.items():
        counter += 1
        if counter % 1000 == 0:
            print(counter)
        gene_seqs = []

        # get protein sequences
        for species in species_list:
            protein_id = orthos[species + ' protein']
            if isinstance(protein_id, float):  # when protein id is nan
                continue
            try:
                gene_seqs.append(fasta_dict[species][protein_id])
            except KeyError:
                nf.write('\t'.join([gene_id, species, protein_id]) + '\n')
        if len(gene_seqs) != 0:
            SeqIO.write(gene_seqs, os.path.join(save_path, "%s.fa" % gene_id), "fasta")
    return


def msa(dir_path, save_path, testn=None):
    fl = [i for i in os.listdir(dir_path) if '.fa' in i]
    if testn:
        fl = fl[:testn]
    nf = open(os.path.join(save_path, 'output.log'), 'w')
    counter = 0
    for f in fl:
        p = os.path.join(dir_path, f)
        sp = os.path.join(save_path, f[:-2] + 'aln')
        if os.path.exists(sp):
            continue
        counter += 1
        if counter % 100 == 0:
            print(counter)
        try:
            subprocess.check_call(['muscle', '-clwstrict', '-in', p, '-out', sp, '-quiet'], stderr=nf)
        except subprocess.CalledProcessError:
            print('Error on', p)
            nf.write(p)
    nf.close()
    return


if __name__ == '__main__':
    # load all peptides fasta
    fasta_dict = get_fasta_dict(fasta_dir='/Volumes/Data/JacJac/genome/')

    # load gene ortho info
    usecols = ['Gene stable ID',
               'Chinese hamster CriGri gene stable ID', 'Chinese hamster CriGri protein or transcript stable ID',
               'Golden Hamster gene stable ID', 'Golden Hamster protein or transcript stable ID',
               'Prairie vole gene stable ID', 'Prairie vole protein or transcript stable ID',
               'Mouse gene stable ID', 'Mouse protein or transcript stable ID',
               'Upper Galilee mountains blind mole rat gene stable ID',
               'Upper Galilee mountains blind mole rat protein or transcript stable ID',
               'Northern American deer mouse gene stable ID',
               'Northern American deer mouse protein or transcript stable ID',
               'Rat gene stable ID', 'Rat protein or transcript stable ID']
    rename_col = ['Jaculus_jaculus gene',
                  'Cricetulus_griseus_crigri gene', 'Cricetulus_griseus_crigri protein',
                  'Mesocricetus_auratus gene', 'Mesocricetus_auratus protein',
                  'Microtus_ochrogaster gene', 'Microtus_ochrogaster protein',
                  'Mus_musculus gene', 'Mus_musculus protein',
                  'Nannospalax_galili gene', 'Nannospalax_galili protein',
                  'Peromyscus_maniculatus_bairdii gene', 'Peromyscus_maniculatus_bairdii protein',
                  'Rattus_norvegicus gene', 'Rattus_norvegicus protein']
    query_col = ['Gene stable ID',
                 'Chinese hamster Query protein or transcript ID',
                 'Golden Hamster Query protein or transcript ID',
                 'Mouse Query protein or transcript ID',
                 'Prairie vole Query protein or transcript ID',
                 'Upper Galilee Query protein or transcript ID',
                 'Northern American deer mouse Query protein or transcript ID',
                 'Rat Query protein or transcript ID']

    #gene_id_dict = get_gene_ids('/Users/hq/Desktop/ortho_total_idonly.tsv', usecols, rename_col,
    #                            query_col=query_col)
    #get_orthogroup_seq(gene_id_dict, fasta_dict, save_path='/Volumes/Data/JacJac/pep/')

    msa('/Volumes/Data/JacJac/pep/', '/Volumes/Data/JacJac/aln/')

