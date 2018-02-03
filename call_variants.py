"""
diff2seq:

block       -   +   +   -   +   +   +   -
true pos    1   1   2   3   4   4   5   6
pos         0   1   2   3   4   5   6   7
ref         A   -   T   K   P   -   C   T
alt         A   -   -   K   C   C   T   T

True pos is the pos count for real aa (not "-"), pos is char count
Identify change block start from pos 0 to the last pos
Block start with first non-match/double-gap pair, end before first match and non-double-gap pair
Then it would be easy to identify var from each block

"""

from Bio import AlignIO
import os
import provean

species_id_dict = {'ENSJJAP': 'Jaculus_jaculus',
                   'ENSCGRP': 'Cricetulus_griseus_crigri',
                   'ENSMAUP': 'Mesocricetus_auratus',
                   'ENSMOCP': 'Microtus_ochrogaster',
                   'ENSMUSP': 'Mus_musculus',
                   'ENSNGAP': 'Nannospalax_galili',
                   'ENSPEMP': 'Peromyscus_maniculatus_bairdii',
                   'ENSRNOP': 'Rattus_norvegicus'}
species_list = ['Cricetulus_griseus_crigri', 'Mesocricetus_auratus', 'Microtus_ochrogaster',
                'Mus_musculus', 'Nannospalax_galili', 'Peromyscus_maniculatus_bairdii', 'Rattus_norvegicus']

write_line = lambda l: '\t'.join(map(str, l)) + '\n'


def check_type(col):
    if col[0] == col[1]:
        if '-' not in col:
            return 'no'
        else:
            return '--'
    elif '-' not in col:
        return 'snv'
    elif '-' in col[0]:
        return 'in'
    else:
        return 'del'


def get_species(gene_id, species_dict=species_id_dict):
    start = gene_id[:7]
    try:
        species = species_dict[start]
    except KeyError:
        species = 'unknown'
    return species


def diff2seq(refseq, altseq, del_char='.'):
    blocks = []
    var_list = []
    true_pos = 0
    block_start = 0
    true_pos_seq = []
    identity_count = 0
    for pos in range(len(refseq)):
        col_type = check_type((refseq[pos], altseq[pos]))
        if refseq[pos] != '-':
            true_pos += 1
        true_pos_seq.append(true_pos)
        if col_type == 'no':
            identity_count += 1
            if block_start == pos:
                block_start = pos + 1
                continue  # no block recording
            else:
                blocks[-1].append(pos)  # add block ending
                block_start = pos + 1
        else:
            if len(blocks) == 0 or len(blocks[-1]) == 2:
                blocks.append([block_start])
            if pos == len(refseq) - 1:
                blocks[-1].append(pos + 1)
    identity = '%.3f' % (identity_count / true_pos)

    for block in blocks:
        block_ref = ''.join([i for i in refseq[block[0]:block[1]] if i != '-'])  # remove '-'
        block_alt = ''.join([i for i in altseq[block[0]:block[1]] if i != '-'])
        block_pos = true_pos_seq[block[0]:block[1]]
        if len(block_ref) == 1 and len(block_alt) == 1:
            block_type = 'sub'
        elif len(block_ref) == 0 and len(block_alt) >= 1:
            block_type = 'ins'
            if block[0] == 0:
                block_ref = str(refseq).lstrip('-')[0]
                block_alt = str(refseq).lstrip('-')[0] + block_alt
                block_type = 'delins'
                # this is actually extension by hgvs, but provean seems to only support this kind expression
            else:
                block_ref = refseq[block[0]-1]
                block_alt = altseq[block[0]-1] + block_alt
        elif len(block_ref) >= 1 and len(block_alt) == 0:
            block_type = 'del'
            block_alt = del_char
        elif len(block_ref) >= 1 and len(block_alt) >= 1:
            block_type = 'delins'
        elif len(block_ref) == 0 and len(block_alt) == 0:
            # both "-", "-" situation
            continue
        else:
            block_type = 'unknown'
            print(block_ref, block_alt, 'unknown type')

        if refseq[block[0]] == '-' and block_type != 'ins':
            # In the case where ref seq block start with "-",
            # e.g. ref "--AA", alt "GCAA",
            # but also prevent to include ins here
            block_pos[0] += 1
        if block_pos[0] == 0:
            block_pos[0] = 1
        block_var = [block_pos[0], block_pos[-1], block_ref, block_alt, block_type]
        var_list.append(block_var)
    return identity, var_list


def call_variants(aln_path, save_dir, key='ENSJJAP', based='key'):
    if os.path.exists(os.path.join(save_dir, os.path.split(aln_path)[1][:-3] + 'var.tsv')):
        return
    align = AlignIO.read(open(aln_path), format='clustal')
    seq_id_list = [i.id.split('.')[0] for i in align._records]
    key_index = -1
    other_index = []
    for i in range(len(seq_id_list)):
        if key in seq_id_list[i]:
            key_index = i
        else:
            other_index.append(i)
    if key_index == -1:
        print(aln_path, key, 'not found.')
        return
    else:
        species_used_list = [get_species(i) for i in seq_id_list]
        del species_used_list[key_index]

    total_var = []
    for other in other_index:
        if based == 'key':
            refseq = align[key_index]
            altseq = align[other]
        elif based == 'ortho':
            refseq = align[other]
            altseq = align[key_index]
        else:
            raise ValueError("Unknown %s based" % based)
        identity, var_list = diff2seq(refseq=refseq, altseq=altseq)
        pre_col = [seq_id_list[key_index], seq_id_list[other], get_species(seq_id_list[other]), identity]
        total_var += [pre_col + i for i in var_list]

    title_ll = ['Jaculus_jaculus protein', 'Orthologous protein', 'Orthologous species', 'Identity',
                'Start', 'End', 'Ref', 'Alt', 'Type']
    nf = open(os.path.join(save_dir, os.path.split(aln_path)[1][:-3] + 'var.tsv'), 'w')
    nf.write(write_line(title_ll))
    for var in total_var:
        nf.write(write_line(var))
    nf.close()


if __name__ == '__main__':
    aln_fl = [os.path.join('/Users/hq/data/jerboa/aln/', i)
              for i in os.listdir('/Users/hq/data/jerboa/aln/') if '.aln' in i]
    counter = 0
    for p in aln_fl:
        counter += 1
        if counter % 100 == 0:
            print(counter)
        try:
            call_variants(p, '/Users/hq/data/jerboa/var/', based='key')
        except ValueError:
            print(p)
    for p in aln_fl:
        counter += 1
        if counter % 100 == 0:
            print(counter)
        try:
            call_variants(p, '/Users/hq/data/jerboa/var_ortho/', based='ortho')
        except ValueError:
            print(p)
        gene_id = os.path.split(p)[1][:-4]
        provean.get_provean_input(gene_id=gene_id, dir_path='/Users/hq/data/jerboa/')
