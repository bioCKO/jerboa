# -*- coding: utf-8 -*-
import dash
import dash_core_components as dcc
import dash_html_components as html
import numpy as np
import pandas as pd
import plotly.figure_factory as ff
from collections import Counter
import dash_table_experiments as dt
from Bio import AlignIO

app = dash.Dash()
# input data, don't modify these inside callback
df = pd.read_table('total_var.filtered.tsv')
protein2gene_series = pd.read_table('jj_gene_protein.txt', sep='\t',
                                    header=0, index_col='Protein stable ID')['Gene stable ID']
provean_score = df['PROVEAN']

color_scheme_Z = {
    # Hydrophobic
    'W': 0.0,
    'L': 0.0,
    'V': 0.0,
    'I': 0.0,
    'M': 0.0,
    'F': 0.0,
    'A': 0.0,
    'C': 0.0,
    # Positive charge
    'K': 0.1,
    'R': 0.1,
    # Negative charge
    'E': 0.2,
    'D': 0.2,
    # Polar
    'N': 0.3,
    'Q': 0.3,
    'S': 0.3,
    'T': 0.3,
    # Cysteines
    'C': 0.4,
    # Glycines
    'G': 0.5,
    # Prolines
    'P': 0.6,
    # Aromatic
    'H': 0.7,
    'Y': 0.7,
    # Other
    '-': 1.0,
}
color_scale = [[0.0, 'rgb(0, 64, 255)'], [0.1, 'rgb(255, 0, 0)'],
               [0.2, 'rgb(191, 0, 255)'], [0.3, 'rgb(0, 255, 0)'],
               [0.4, 'rgb(255, 0, 255)'], [0.5, 'rgb(255, 128, 0)'],
               [0.6, 'rgb(255, 255, 0)'], [0.7, 'rgb(0, 191, 255)'],
               [1.0, 'rgb(255, 255, 255)']]
species_id_dict = {'ENSJJAP': 'Jerboa',
                   'ENSCGRP': 'Chinese hamster',
                   'ENSMAUP': 'Golden Hamster',
                   'ENSMOCP': 'Prairie vole',
                   'ENSMUSP': 'Mouse',
                   'ENSNGAP': 'Blind mole rat',
                   'ENSPEMP': 'Deer mouse',
                   'ENSRNOP': 'Rat'}


def trans_ppid_species(ppl, id_dict=species_id_dict):
    return [id_dict[i[:7]] for i in ppl]


def aln_color_z(aln, color_code=color_scheme_Z):
    shape = aln.shape
    coloring_count_threshold = int(min(shape[0]-1, shape[0]*0.8))
    color_cols = []
    for i in range(shape[1]):
        cur_col = aln[:, i]
        cur_col_count = Counter(cur_col)
        if len(sorted(cur_col_count.elements())) == 1: # all the same
            coloring = [cur_col_count.elements()]
        else:
            coloring = []
            for j in cur_col_count.most_common(3):
                if j[1] >= coloring_count_threshold:
                    coloring.append(j[0])
        cur_color = []
        for aa in cur_col:
            if aa in coloring and aa != '-':
                cur_color.append(color_code[aa])
            else:
                cur_color.append(color_code['-'])
        color_cols.append(cur_color)
    z = np.array(color_cols).T.tolist()
    return z


def true_pos2msa_pos(msa_seq, query_true_pos):
    pos_count = 0
    query_pos = 0
    true_pos_list = []
    for i in range(len(msa_seq)):
        if msa_seq[i] != '-':
            pos_count += 1
        if pos_count == query_true_pos:
            query_pos = i
        true_pos_list.append(str(pos_count))
    return query_pos, true_pos_list


app.layout = html.Div(children=[
    # Hidden div inside the app that stores the intermediate value
    html.Div(id='intermediate-value', style={'display': 'none'}),

    html.H1(children='Jerboa Protein Variants Analysis'),
    dcc.Graph(
        id='PROVEAN-score-distribution',
    ),
    dcc.Slider(
        id='threshold-slider',
        min=provean_score.min(),
        max=provean_score.max() - 0.1,
        value=provean_score.quantile(.01),
        step=0.05,
        updatemode='mouseup'
    ),
    html.Div(children=[
        dcc.Graph(id='pro-var-count'),
        html.H3(id='stat-str'),
    ]),
    html.Div(children=[
        html.H4('Variants Table'),
        dt.DataTable(
            rows=[{}],
            columns=df.columns,
            row_selectable=True,
            filterable=True,
            sortable=True,
            selected_row_indices=[0],
            editable=False,
            id='var-table'
        ),
    ]),
    html.Div(id='var-info', children=[
        html.Div(id='var-str'),
        dcc.Graph(id='var-msa-matrix')
    ])
])


@app.callback(
    dash.dependencies.Output('PROVEAN-score-distribution', 'figure'),
    [dash.dependencies.Input('threshold-slider', 'value')]
)
def update_distplot(threshold):
    hist_data = [provean_score]
    group_labels = ['']
    fig = ff.create_distplot(hist_data, group_labels, bin_size=0.2,
                             show_hist=True, show_curve=True, show_rug=False)
    fig.layout['title'] = 'Overall PROVEAN Score Distribution'
    fig.layout['shapes'] = [
        # Line reference to the axes
        {
            'type': 'rect',
            'xref': 'x',
            'yref': 'paper',
            'x0': threshold,
            'y0': 0,
            'x1': provean_score.max(),
            'y1': 1,
            'fillcolor': '#FFFFFF',
            'opacity': 0.5,
            'line': {
                'width': 0,
            }
        },
        {
            'type': 'line',
            'xref': 'x',
            'yref': 'y',
            'x0': threshold,
            'y0': 0,
            'x1': threshold,
            'y1': 0.5,
            'line': {
                'color': 'rgb(55, 128, 191)',
                'width': 3,
            },
        },
    ]
    return fig


@app.callback(
    dash.dependencies.Output('intermediate-value', 'children'),
    [dash.dependencies.Input('threshold-slider', 'value')]
)
def filter_var_by_provean_threshold(threshold):
    sub_df = df[provean_score < threshold]
    return sub_df.to_json(date_format='iso', orient='split')


@app.callback(
    dash.dependencies.Output('pro-var-count', 'figure'),
    [dash.dependencies.Input('intermediate-value', 'children')]
)
def update_var_count_bar(data):
    sub_df = pd.read_json(data, orient='split')
    jj_protein_var_counts = sub_df['Jaculus_jaculus protein'].value_counts()
    protein_var_stat = pd.cut(jj_protein_var_counts,
                              bins=[0, 1, 2, 3, 5, 10, 25, 1000], right=True,
                              labels=['Exact 1', 'Exact 2', 'Exact 3', '4-5', '6-10', '11-25', '>25']).value_counts()
    fig = {
        'data': [
            {'x': protein_var_stat.index.tolist(),
             'y': protein_var_stat.tolist(),
             'type': 'bar', 'name': 'Count'},
        ],
        'layout': {
            'title': 'Number of Variants Remained in Each Gene'

        }
    }
    return fig


@app.callback(
    dash.dependencies.Output('stat-str', 'children'),
    [dash.dependencies.Input('intermediate-value', 'children')]
)
def update_var_gene_num(data):
    sub_df = pd.read_json(data, orient='split')
    jj_protein_var_counts = sub_df['Jaculus_jaculus protein'].value_counts()
    total_var = jj_protein_var_counts.sum()
    total_gene = len(jj_protein_var_counts)
    return 'Total Variants Remained %d, Total Genes Remained %d' % (total_var, total_gene)


@app.callback(
    dash.dependencies.Output('var-table', 'rows'),
    [dash.dependencies.Input('intermediate-value', 'children')]
)
def update_var_table(data):
    sub_df = pd.read_json(data, orient='split')
    return sub_df.to_dict('records')


@app.callback(
    dash.dependencies.Output('var-msa-matrix', 'figure'),
    [dash.dependencies.Input('var-table', 'selected_row_indices')],
    [dash.dependencies.State('intermediate-value', 'children')]
)
def update_single_var_display(selected_row_indices, data):
    sub_df = pd.read_json(data, orient='split')
    if len(selected_row_indices) == 0:
        return 'Nothing selected.'
    var = sub_df.iloc[selected_row_indices[0]]  # the Series of last selected var
    jj_protein_id = var['Jaculus_jaculus protein']
    jj_true_pos = var['Start']
    jj_gene_id = protein2gene_series[jj_protein_id]
    jj_aln = AlignIO.read(open('aln/%s.aln' % jj_gene_id), format='clustal')
    seq_id_list = [i.id.split('.')[0] for i in jj_aln._records]
    jj_aln_idx = seq_id_list.index(jj_protein_id)
    jj_length = len(jj_aln[jj_aln_idx])
    jj_pos, true_pos_list = true_pos2msa_pos(jj_aln[jj_aln_idx], jj_true_pos)  # change this to some func
    region_start = max(0, jj_pos-15)
    region_end = min(jj_length, jj_pos + 15)
    plot_jj_aln = jj_aln[:, region_start:region_end]
    x = true_pos_list[region_start:region_end]
    y = trans_ppid_species([r.id for r in plot_jj_aln])

    symbol = np.array(plot_jj_aln).tolist()
    z = aln_color_z(np.array(plot_jj_aln))
    fig = ff.create_annotated_heatmap(z, x=None, y=y, annotation_text=symbol, colorscale=color_scale,
                                      font_colors=['black'], hoverinfo="x+y")
    fig.layout.xaxis.tickmode = 'array'
    fig.layout.xaxis.tickvals = list(range(len(x)))
    fig.layout.xaxis.ticktext = x
    fig.layout.margin.l = 120
    shape_x_add = min(14, var['End']-var['Start'])
    fig.layout.shapes = [
        {   # shape for the variants
            'type': 'rect',
            'x0': 14.5 if jj_pos >= 15 else jj_pos,
            'y0': -1,
            'x1': 15.5 + shape_x_add if jj_pos >= 15 else jj_pos,
            'y1': len(seq_id_list)+0.5,
            'line': {
                'color': 'rgba(128, 128, 128, 0.8)',
                'width': 2,
            },
            'fillcolor': 'rgba(128, 128, 128, 0.4)',
        },
        {  # shape for jerboa sequence
            'type': 'rect',
            'x0': -0.5,
            'y0': jj_aln_idx - 0.5,
            'x1': region_end-region_start-0.5,
            'y1': jj_aln_idx + 0.5,
            'line': {
                'color': 'rgba(0, 191, 255, 1)',
                'width': 2,
            },
            'fillcolor': 'rgba(0, 191, 255, 0.5)',
        }
    ]
    return fig


@app.callback(
    dash.dependencies.Output('var-str', 'children'),
    [dash.dependencies.Input('var-table', 'selected_row_indices')],
    [dash.dependencies.State('intermediate-value', 'children')]
)
def update_var_char(selected_row_indices, data):
    sub_df = pd.read_json(data, orient='split')
    if len(selected_row_indices) == 0:
        return 'Nothing selected.'
    var = sub_df.iloc[selected_row_indices[0]]  # the Series of last selected var
    value = [html.H3('Selected Variant Information'),
             html.P('Jerboa Protein ID: {jjp}'.format(jjp=var['Jaculus_jaculus protein'])),
             html.P('Orthologous Protein ID: {othop}'.format(othop=var['Orthologous_protein'])),
             html.P('Orthologous Protein to Jerboa Protein: {idt}'.format(idt=var['Identity'])),
             html.P('Variant Impact Region: {start}-{end}'.format(start=var['Start'], end=var['End'])),
             html.P('Ref: {ref}; Alt: {alt}'.format(ref=var['Ref'], alt=var['Alt'])),
             html.P('PROVEAN score: {provean}'.format(provean='%.3f' % var['PROVEAN']))]
    return value


if __name__ == '__main__':
    app.run_server(debug=True, port=5063)
