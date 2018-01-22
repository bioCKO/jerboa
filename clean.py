import hashlib
import pandas as pd
import os
import warnings
from multiprocessing import Pool

warnings.filterwarnings('ignore')

path = '/Volumes/Data/database/ncbi/nr'


def check_md5(path):
    fl = os.listdir(path)
    nf = open(path + '.check.txt', 'w')

    for f in fl:
        if 'nr' in f and 'md5' not in f:
            print(f)
            md5 = hashlib.md5(open(os.path.join(path, f), 'rb').read()).hexdigest()
            real_md5 = open(os.path.join(path, f) + '.md5').read().split(' ')[0]
            if md5 == real_md5:
                nf.write('True' + '  ' + f + '\n')
                print('True' + '  ' + f + '\n')
            else:
                nf.write('False' + ' ' + f + '\n')
                print('False' + '  ' + f + '\n')
    nf.close()
    return


def select_ortho(path):
    """
    merge lines when gene have multiple transcripts isoform
    select lines when gene have multiple orthologous genes, chose the one with max identity% to query gene.
    in the result table, each line represents a gene, and its unique ortholougous gene (if exist) in other species.
    the number of lines equals to the number of genes in the query species.
    :param path:
    :return:
    """
    print(path)
    univer_title = ['Gene stable ID', 'Transcript stable ID', 'Protein stable ID',
                    'Chromosome/scaffold name', 'Gene start (bp)', 'Gene end (bp)', 'Strand',
                    'Karyotype band', 'Gene name', 'Source of gene name', 'Transcript count',
                    'Gene % GC content', 'Gene description', 'Target gene stable ID', 'Target gene name',
                    'Target protein or transcript stable ID', 'Target chromosome/scaffold name',
                    'Target chromosome/scaffold start (bp)', 'Target chromosome/scaffold end (bp)',
                    'Query protein or transcript ID', 'Last common ancestor with Target',
                    'Target homology type', '%id. target gene identical to query gene',
                    '%id. query gene identical to target gene', 'Target Gene-order conservation score',
                    'Target Whole-genome alignment coverage', 'dN with Target', 'dS with Target',
                    'Target orthology confidence [0 low, 1 high]']
    df = pd.read_table(path)
    df.set_index('Gene stable ID', inplace=True, drop=False)
    unique_geneid = df.index.unique()
    original_title = df.columns
    df.columns = univer_title
    merge_info = lambda i: '|'.join(i.unique().tolist())
    total_lines = []
    counter = 0

    for geneid in unique_geneid:
        counter += 1
        if counter % 2500 == 0:
            print(counter)
        gene_info = df.loc[geneid, :]
        if isinstance(gene_info, pd.Series):
            total_lines.append(gene_info.tolist())
        else:
            gene_info_line = gene_info.iloc[0]
            gene_info_line['Transcript stable ID'] = merge_info(gene_info['Transcript stable ID'])
            try:
                gene_info_line['Protein stable ID'] = merge_info(gene_info['Protein stable ID'])
            except TypeError:
                pass
            if gene_info_line['Target homology type'] != 'ortholog_one2one' and not \
                    pd.isnull(gene_info['Target gene stable ID']).all():
                try:
                    gene_info_line_target = gene_info[gene_info['%id. target gene identical to query gene'] \
                                                      == gene_info['%id. target gene identical to query gene'] \
                                                          .max()].iloc[0].tolist()
                except IndexError:
                    print(gene_info[gene_info['%id. target gene identical to query gene'] \
                                    == gene_info['%id. target gene identical to query gene'] \
                          .max()].shape)
                    print(gene_info)
                    raise
                new_line = gene_info_line.tolist()[:13] + gene_info_line_target[13:]
            else:
                new_line = gene_info_line.tolist()
            total_lines.append(new_line)
    pd.DataFrame(total_lines, columns=original_title).to_csv(path[:-3] + 'select.tsv', sep='\t', index=False)
    print(path)
    return


if __name__ == '__main__':
    ortho_table_dir = '/Volumes/Data/kim_lab/ortho/'
    fl = os.listdir(ortho_table_dir)[1:]

    pool = Pool(8)
    for f in fl:
        print(f)
        pool.apply_async(select_ortho, (ortho_table_dir + f,))
    pool.close()
    pool.join()




