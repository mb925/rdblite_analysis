from functools import reduce
from scipy.stats import pearsonr
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix
import pandas as pd
import json
import config as cfg
from pymongo import MongoClient


# create input.txt file from entries
def create_input_txt():
    f = open(cfg.data['data'] + '/entries.json')
    data = json.load(f)
    pdb_list = []
    for i in data:
        if i['repeatsdb_id'].startswith('4acqu') or i['repeatsdb_id'] == '5b7oA': # excluding obsoleted PDBs
            continue
        if i['reviewed'] == True:
            print(i['reviewed'])
            pdb_list.append(i['repeatsdb_id'])

    print(pdb_list)
    pdb_list = set(pdb_list)

    with open(cfg.data['data'] + '/input.txt', 'a') as the_file:
        the_file.write('pdb_id' + ' ' + 'pdb_chain' + '\n')
        for el in pdb_list:
            the_file.write(el[0:4] + ' ' + el[4] + '\n')

def create_input_negatives():
    df = pd.read_csv(cfg.data['data'] + '/non_trp_pdbchains.tsv', sep='\t')
    df = df.drop(columns=["UniProtID"])
    df['pdb_id'], df['pdb_chain'] = df['pdb_chain'].str[4:], df['pdb_chain'].str[:4]
    df.to_csv(cfg.data['data'] + "/input_negatives.txt", index=False, sep=' ')
    print(df)
def table_class_topology():

    f = open(cfg.data['data'] + '/entries.json')
    pdbs_list = pd.read_csv(cfg.data['data'] + '/input.txt', sep=' ')
    data = json.load(f)
    data = pd.DataFrame.from_records(data)
    data = data.loc[data['class'] != ""]
    data = pd.merge(data, pdbs_list, on=['pdb_id', 'pdb_chain'], how='right')

    # pdb chains
    pdbs_top = data.groupby(['repeatsdb_id','class_topology'])['class_topology'].size().reset_index(name='# of PDB chains')
    pdbs_top = pdbs_top.groupby(['class_topology']).size().reset_index(name='# of PDB chains')
    pdbs_top = pdbs_top.append(pdbs_top.sum(numeric_only=True), ignore_index=True).fillna('Total')
    pdbs_cl = data.groupby(['repeatsdb_id','class'])['class'].size().reset_index(name='# of PDB chains')
    pdbs_cl = pdbs_cl.groupby(['class']).size().reset_index(name='# of PDB chains')
    pdbs_cl = pdbs_cl.append(pdbs_cl.sum(numeric_only=True), ignore_index=True).fillna('Total')

    # regions
    regions_top = data.groupby(['region_id','class_topology'])['class_topology'].size().reset_index(name='regions_count')
    regions_top = regions_top.groupby(['class_topology']).size().reset_index(name='regions_count')
    regions_top = regions_top.append(regions_top.sum(numeric_only=True), ignore_index=True).fillna('Total')
    regions_cl = data.groupby(['region_id','class'])['class'].size().reset_index(name='regions_count')
    regions_cl = regions_cl.groupby(['class']).size().reset_index(name='regions_count')
    regions_cl = regions_cl.append(regions_cl.sum(numeric_only=True), ignore_index=True).fillna('Total')

    # units
    units = data.loc[data['type'] == 'unit']
    units_top = units.groupby(['class_topology'])['class_topology'].size().reset_index(name='units_count')
    units_top = units_top.append(units_top.sum(numeric_only=True), ignore_index=True).fillna('Total')
    units_cl = units.groupby(['class'])['class'].size().reset_index(name='units_count')
    units_cl = units_cl.append(units_cl.sum(numeric_only=True), ignore_index=True).fillna('Total')

    # avg_units_len
    avg_units_top = units.groupby(['class_topology'])['region_average_unit_length'].mean().round(2).reset_index(name='avg_units')
    avg_units_top = avg_units_top.append(avg_units_top.mean().round(2), ignore_index=True).fillna('Total')
    avg_units_cl = units.groupby(['class'])['region_average_unit_length'].mean().round(2).reset_index(name='avg_units')
    avg_units_cl = avg_units_cl.append({'class': 'Total', 'avg_units': avg_units_cl['avg_units'].mean().round(2)}, ignore_index=True)


    # std
    std_units_top = units.groupby(['class_topology'])['region_average_unit_length'].std().round(2).reset_index(
        name='std_units')
    std_units_top = std_units_top.append({'class_topology': 'Total', 'std_units': std_units_top['std_units'].mean().round(2)}, ignore_index=True)

    std_units_cl = units.groupby(['class'])['region_average_unit_length'].std().round(2).reset_index(
        name='std_units')
    std_units_cl = std_units_cl.append({'class': 'Total', 'std_units': std_units_top['std_units'].mean().round(2)}, ignore_index=True).fillna('Total') # todo append total

    # uniprot
    uniprot = pd.read_csv(cfg.data['data'] + '/pdb_chain_uniprot.tsv',  sep='\t', header=1, index_col='SP_PRIMARY').reset_index()
    uniprot['repeatsdb_id']=uniprot['PDB'] + uniprot['CHAIN']
    up = uniprot[['SP_PRIMARY', 'repeatsdb_id']]

    uniprot_pdb = data.merge(up, how='left', on='repeatsdb_id')
    uniprot_top = uniprot_pdb.groupby(['SP_PRIMARY','class_topology'])['class_topology'].size().reset_index(name='uniprot_counts')
    uniprot_top = uniprot_top.groupby(['class_topology'])['class_topology'].size().reset_index(name='uniprot_counts')
    uniprot_top = uniprot_top.append(uniprot_top.sum(numeric_only=True), ignore_index=True).fillna('Total')
    uniprot_cl = uniprot_pdb.groupby(['SP_PRIMARY', 'class'])['class'].size().reset_index(name='uniprot_counts')
    uniprot_cl = uniprot_cl.groupby(['class'])['class'].size().reset_index(name='uniprot_counts')
    uniprot_cl = uniprot_cl.append(uniprot_cl.sum(numeric_only=True), ignore_index=True).fillna('Total')

    class_table = reduce(lambda x, y: pd.merge(x, y, on='class', how='outer'), [pdbs_cl, regions_cl, units_cl, avg_units_cl, std_units_cl, uniprot_cl])
    topology_table = reduce(lambda x, y: pd.merge(x, y, on='class_topology', how='outer'), [pdbs_top, regions_top, units_top, avg_units_top, std_units_top, uniprot_top])

    topology_table.to_csv(cfg.data['data'] + "/topologies_count.csv")
    class_table.to_csv(cfg.data['data'] + "/classes_count.csv")
    print(pdbs_top)

def residues_count():
    data = pd.read_csv(cfg.data['data'] + '/binary_pdb.csv', sep=',')
    repeated = len(data.loc[data['CURATED'] == 1])
    not_repeated = len(data.loc[data['CURATED'] == 0])
    fig = go.Figure()
    fig.add_trace(go.Bar(
        y=['RepeatsDB-lite 2 TRP', 'RepeatsDB-lite 1 TRP', 'Curated TRP'],
        x=[len(data.loc[data['RDB2'] == 1]) , len(data.loc[data['RDB1'] == 1]), repeated],
        name='Repeated residued',
        orientation='h',
        marker=dict(
            color='rgba(246, 78, 139, 0.6)',
            line=dict(color='rgba(246, 78, 139, 1.0)', width=3)
        )
    ))
    fig.add_trace(go.Bar(
        y=['RepeatsDB-lite 2 TRP', 'RepeatsDB-lite 1 TRP', 'Curated TRP'],
        x=[len(data.loc[data['RDB2'] == 0]), len(data.loc[data['RDB1'] == 0]), not_repeated],
        name='Not repeated',
        orientation='h',
        marker=dict(
            color='rgba(58, 71, 80, 0.6)',
            line=dict(color='rgba(58, 71, 80, 1.0)', width=3)
        )
    ))

    fig.update_layout(barmode='stack', height=300)
    fig.write_image(cfg.data['plots'] + "/residues_count.png")

def parse_curated_to_dict():
    print('entries dict')
    dict_curated = {}
    f = open(cfg.data['data'] + '/entries.json')
    data = json.load(f)
    for i in data:
        if i['repeatsdb_id'].startswith('4acqu') or i['repeatsdb_id'] == '5b7oA': # excluding obsoleted PDBs
            continue
        if i['reviewed'] == True:
            if i["repeatsdb_id"] in dict_curated:
                dict_curated[i["repeatsdb_id"]].append([i['start'], i['end']])
            if i["repeatsdb_id"] not in dict_curated:
                dict_curated[i["repeatsdb_id"]] = []
                dict_curated[i["repeatsdb_id"]].append([i['start'], i['end']])
    print(dict_curated['4acqA'])
    return dict_curated


def parse_regions_to_dict():
    print('regions dict')
    dict_regions = {}
    f = open(cfg.data['data'] + '/entries.json')
    data = json.load(f)
    for i in data:
        if i['reviewed'] == True:
            if i["repeatsdb_id"] in dict_regions:
                st = i["region_id"].split('_')[1]
                end = i["region_id"].split('_')[2]
                dict_regions[i["repeatsdb_id"]].append([st, end])
            if i["repeatsdb_id"] not in dict_regions:
                st = i["region_id"].split('_')[0]
                end = i["region_id"].split('_')[1]
                dict_regions[i["repeatsdb_id"]] = []
                dict_regions[i["repeatsdb_id"]].append([st, end])

    return dict_regions


def parse_predicted2_to_dict(f):
    print('predicted dict')
    dict_predicted = {}
    data = pd.read_csv(f, sep=',')

    for i in data.iterrows():
        if i[1][2] == 'no regions':
            continue
        units = json.loads(i[1][2])
        if i[1][0] in dict_predicted:
            dict_predicted[i[1][0]].extend(units)
        if i[1][0] not in dict_predicted:
            dict_predicted[i[1][0]] = []
            dict_predicted[i[1][0]].extend(units)
    return dict_predicted


def parse_predicted1_to_dict(f):
    print('predicted dict')
    dict_predicted = {}
    df = pd.read_csv(f, sep='\t', names=["PDB", "TYPE", "START", "END"])
    df_intervals = df.groupby('PDB').apply(lambda x: build_intervals(x))
    for row in df_intervals:
        dict_predicted[row[0]] = row[1]
    return dict_predicted


def build_intervals(x):
    intervals = []
    units = x.loc[x['TYPE'] == 'UNIT']
    # print(x)
    pdb = units.iloc[0]['PDB']
    pdb = pdb[3:].split('_')[0] + pdb[3:].split('_')[1]
    for row in units.iterrows():
        intervals.append([str(row[1]['START']), str(row[1]['END'])])
    return [pdb, intervals]


def find_seqres(db, pdb_id, pdb_chain):

    print('Current pdb is %s' % str(pdb_id))
    print('Current pdb is %s' % str(pdb_chain))
    l = list(db.seqres.find({'pdb_id': pdb_id, 'pdb_chain': pdb_chain}))
    return l


def dict_to_binary(d_cur, d_pred1, d_pred2, d_reg):
    # Connect to mongo (BioDBs seqres)
    client = MongoClient(cfg.db_host, cfg.db_port)
    db = client.biodbs


    df = pd.DataFrame()

    for pdb in d_cur:

        l = find_seqres(db, pdb[0:4], pdb[-1])
        df_residues = pd.DataFrame.from_records(l)
        if df_residues.empty:
            print('empty ', pdb)
            continue
        residues = df_residues.sort_values(by=['seqres_index'], ascending=True)[
            ['pdb_id', 'pdb_chain', 'seqres_index', 'pdb_residue_id']]
        residues['REGION'] = '0'
        residues['CURATED'] = 0
        residues['RDB1'] = 0
        residues['RDB2'] = 0
        for interval in d_cur[pdb]:
            fill_repeated_interval(pdb, residues, interval, 'CURATED')  # fill with 1 in the repeated interval
        if pdb in d_pred1:
            for interval1 in d_pred1[pdb]:
                fill_repeated_interval(pdb, residues, interval1, 'RDB1')
        if pdb in d_pred2:
            for interval2 in d_pred2[pdb]:
                fill_repeated_interval(pdb, residues, interval2, 'RDB2')
        if pdb in d_reg:
            for interval3 in d_reg[pdb]:
                fill_repeated_interval(pdb, residues, interval3, 'REGION')

        df = pd.concat([df, residues])
    df["PDB"] = df['pdb_id'] + df["pdb_chain"]
    df.to_csv(cfg.data['data'] + '/binary_pdb.csv', index=False)


def fill_repeated_interval(pdb, residues, interval, type):
    start_seqres = residues.loc[residues['pdb_residue_id'] == str(interval[0])]
    end_seqres = residues.loc[residues['pdb_residue_id'] == str(interval[1])]
    if not start_seqres.empty:
        start_seqres = start_seqres['seqres_index'].iloc[0]
        if not end_seqres.empty:
            end_seqres = end_seqres['seqres_index'].iloc[0]
            if type == 'REGION':
                residues.loc[(residues['seqres_index'] >= start_seqres) &
                             (residues['seqres_index'] <= end_seqres), type] = pdb + '_' + str(interval[0]) + '_' + str(
                    interval[1])
            else:
                residues.loc[(residues['seqres_index'] >= start_seqres) &
                             (residues['seqres_index'] <= end_seqres), type] = 1
    return residues


def calc_hits(x):
    tp1, tn1, fp1, fn1 = 0, 0, 0, 0
    tp2, tn2, fp2, fn2 = 0, 0, 0, 0
    for row in x.iterrows():
        if row[1]['CURATED'] == 1 and row[1]['RDB1'] == 1:
            tp1 += 1
        if row[1]['CURATED'] == 1 and row[1]['RDB1'] == 0:
            fn1 += 1
        if row[1]['CURATED'] == 0 and row[1]['RDB1'] == 1:
            fp1 += 1
        if row[1]['CURATED'] == 0 and row[1]['RDB1'] == 0:
            tn1 += 1
        if row[1]['CURATED'] == 1 and row[1]['RDB2'] == 1:
            tp2 += 1
        if row[1]['CURATED'] == 1 and row[1]['RDB2'] == 0:
            fn2 += 1
        if row[1]['CURATED'] == 0 and row[1]['RDB2'] == 1:
            fp2 += 1
        if row[1]['CURATED'] == 0 and row[1]['RDB2'] == 0:
            tn2 += 1
        print(row[1]['PDB'])

    return [tp1, tn1, fp1, fn1, tp2, tn2, fp2, fn2]


def calc_metrics(x):
    hits = calc_hits(x)
    tp1, tn1, fp1, fn1 = hits[0], hits[1], hits[2], hits[3]
    tp2, tn2, fp2, fn2 = hits[4], hits[5], hits[6], hits[7]

    # balanced accuracy
    se1_total = tp1 + fn1
    se2_total = tp2 + fn2
    sp1_total = fp1 + tn1
    sp2_total = fp2 + tn2
    if se1_total == 0 or sp1_total == 0:
        accuracy1 = 0
    else:
        accuracy1 = ((tp1 / se1_total) + (tn1 / sp1_total)) / 2

    if se2_total == 0 or sp2_total == 0:
        accuracy2 = 0
    else:
        accuracy2 = ((tp2 / se2_total) + (tn2 / sp2_total)) / 2


    p1_total = tp1 + fp1
    p2_total = tp2 + fp2
    r1_total = tp2 + fn2
    r2_total = tp2 + fn2

    if p1_total == 0:
        precision1 = 0
    else:
        precision1 = round(tp1 / p1_total, 2)
    if p2_total == 0:
        precision2 = 0
    else:
        precision2 = round(tp2 / p2_total, 2)
    if r1_total == 0:
        recall1 = 0
    else:
        recall1 = round(tp1 / r1_total, 2)
    if r2_total == 0:
        recall2 = 0
    else:
        recall2 = round(tp2 / r2_total, 2)

    if precision1 == 0 and recall1 == 0:
        f_score1 = 0
    else:
        f_score1 = 2 * ((precision1 * recall1)/(precision1+recall1))
    if precision2 == 0 and recall2 == 0:
        f_score2 = 0
    else:
        f_score2 = 2 * ((precision2 * recall2)/(precision2+recall2))
    print(accuracy1)
    return [accuracy1, accuracy2, precision1, precision2, recall1, recall2, f_score1, f_score2]


def calc_fractions():
    data = pd.read_csv(cfg.data['data'] + '/binary_pdb.csv', sep=',')
    metrics = calc_metrics(data)
    analysis = [[metrics[0], metrics[1], metrics[2], metrics[3], metrics[4], metrics[5], metrics[6], metrics[7]]]
    df = pd.DataFrame(analysis,
                      columns=['accuracy rdbl1', 'accuracy rdbl2', 'precision rdbl1', 'precision rdbl2', 'recall rdbl1',
                               'recall rdbl2', 'mcc1', 'mcc2'])
    df.to_csv(cfg.data['data'] + '/analysis.csv', index=False)


def calc_fractions_pdbs():
    df = pd.read_csv(cfg.data['data'] + '/binary_pdbs_classes.csv', sep=',')

    df = df[df['classes'].notna() & df['topologies'].notna()]
    df = df[~df["topologies"].astype(str).str.contains(';')]

    metrics = df.groupby(['PDB','topologies']).apply(lambda x: calc_metrics(x)).reset_index()
    metrics.rename(columns={0: 'METRICS'}, inplace=True)
    metrics.assign(METRICS=metrics.METRICS.str.split(",")).explode('METRICS')
    metrics2 = pd.DataFrame(metrics.METRICS.values.tolist()).add_prefix('m_')
    metrics2 = metrics2.join(metrics[['PDB', 'topologies']])
    metrics_final = metrics2.rename(
        columns={'m_0': 'a1', 'm_1': 'a2', 'm_2': 'p1', 'm_3': 'p2', 'm_4': 'r1', 'm_5': 'r2', 'm_6': 'f1',
                 'm_7': 'f2'})

    metrics_final.to_csv(cfg.data['data'] + '/analysis_pdbs.csv', index=False)

def calc_fractions_dataset():
    df = pd.read_csv(cfg.data['data'] + '/binary_pdbs_classes.csv', sep=',')

    # remove duplicates
    df = df[df['classes'].notna() & df['topologies'].notna()]
    df = df[~df["topologies"].astype(str).str.contains(';')]

    metrics = df.groupby('topologies').apply(lambda x: calc_metrics(x)).reset_index()
    metrics.rename(columns={0: 'METRICS'}, inplace=True)
    metrics.assign(METRICS=metrics.METRICS.str.split(",")).explode('METRICS')
    metrics2 = pd.DataFrame(metrics.METRICS.values.tolist()).add_prefix('m_')
    metrics2 = metrics2.join(metrics['topologies'])
    metrics_final = metrics2.rename(
        columns={'m_0': 'a1', 'm_1': 'a2', 'm_2': 'p1', 'm_3': 'p2', 'm_4': 'r1', 'm_5': 'r2', 'm_6': 'f1',
                 'm_7': 'f2'})

    metrics_final.to_csv(cfg.data['data'] + '/analysis_classes.csv', index=False)
def plot_f_score():
    df = pd.read_csv(cfg.data['data'] + '/analysis_pdbs.csv', sep=',', dtype={'topologies': str})
    df_mcc = df[['f1', 'f2','topologies']]
    df_mcc = df_mcc[df_mcc['topologies'].notna()]
    df_mcc = df_mcc[~df_mcc["topologies"].str.contains(';')]  # only 1 class
    fig = px.scatter(df_mcc, x='f1', y='f2', color='topologies', title="F-score",
                     height=600)
    fig.update_traces(marker={'size': 3})
    fig.update_layout(xaxis_title="rdbl1", yaxis_title="rdbl2")
    corr, _ = pearsonr(df_mcc['f1'], df_mcc['f2'])
    print('Pearsons correlation: %.3f' % corr)

    fig.write_image(cfg.data['plots'] + "/f-score_topologies.png")


def plot_table():
    df = pd.read_csv(cfg.data['data'] + '/analysis.csv', sep=',')
    fig = go.Figure(data=[go.Table(
        header=dict(values=['predictor', 'accuracy', 'precision', 'recall', 'mcc']),
        cells=dict(values=[['rdblite1', 'rdblite2'], [df['accuracy rdbl1'], df['accuracy rdbl2']],
                           [df['precision rdbl1'], df['precision rdbl2']], [df['recall rdbl1'], df['recall rdbl2']],
                           [df['mcc1'], df['mcc2']]]))
    ])
    fig.write_image(cfg.data['plots'] + "/accuracy_precision_recall.png")


def plot_matrix():
    # Calculate the confusion matrix

    df = pd.read_csv(cfg.data['data'] + '/binary_pdb.csv', sep=',')
    conf_matrix = confusion_matrix(y_true=df['CURATED'], y_pred=df['RDB2'])
    # Print the confusion matrix using Matplotlib
    #
    fig, ax = plt.subplots(figsize=(5, 5))
    ax.matshow(conf_matrix, cmap=plt.cm.Oranges, alpha=0.3)
    for i in range(conf_matrix.shape[0]):
        for j in range(conf_matrix.shape[1]):
            ax.text(x=j, y=i, s=conf_matrix[i, j], va='center', ha='center', size='xx-large')

    plt.xlabel('Predictions', fontsize=18)
    plt.ylabel('Actuals', fontsize=18)
    plt.title('Confusion Matrix - RDBL2', fontsize=18)
    plt.savefig(cfg.data['plots'] + "/matrix_rdbl2.png")
    print(conf_matrix)


def extract_particular_cases():
    df = pd.read_csv(cfg.data['data'] + '/analysis_pdbs.csv', sep=',')
    df1 = df.loc[(df['f2'] == 0)]['PDB']
    df1.to_csv(cfg.data['data'] + '/particular_cases/rdb1_f1_0.csv', index=False)

    # print(df)


def curated_to_dict_classes(level):
    print('entries dict')
    dict_curated = {}
    f = open(cfg.data['data'] + '/entries.json')
    data = json.load(f)
    for i in data:
        if i["repeatsdb_id"] in dict_curated:
            if i[level] not in dict_curated[i["repeatsdb_id"]]:
                dict_curated[i["repeatsdb_id"]].append(i[level])
        if i["repeatsdb_id"] not in dict_curated:
            dict_curated[i["repeatsdb_id"]] = []
            if i[level] not in dict_curated[i["repeatsdb_id"]]:
                dict_curated[i["repeatsdb_id"]].append(i[level])
    return dict_curated


def add_class_to_binary(dict_classes, dict_classes2):
    df = pd.read_csv(cfg.data['data'] + '/binary_pdb.csv', sep=',')
    classes = []
    topologies = []
    for row in df.iterrows():
        print(row[1]['PDB'])
        classes.append(';'.join(map(str, dict_classes[row[1]['PDB']])))
        topologies.append(';'.join(map(str, dict_classes2[row[1]['PDB']])))
    df['classes'] = classes
    df['topologies'] = topologies
    df.to_csv(cfg.data['data'] + '/binary_pdbs_classes.csv', index=False)


def plot_metrics_dataset():
    df = pd.read_csv(cfg.data['data'] + '/analysis_classes.csv', sep=',', dtype={'topologies': str})
    x = ['2.1', '2.2', '3.1', '3.2', '3.3', '3.4',
         '3.6', '4.1', '4.2', '4.3', '4.4', '4.5',
         '4.6', '4.7', '4.8', '4.9', '4.10',
         '5.1', '5.2', '5.3', '5.4', '5.5']
    df = df.sort_values(by=['topologies'], ascending=True)
    # RDBL 1
    fig = go.Figure()


    fig.add_trace(go.Bar(
        y=df['a1'].values.tolist(),
        x=x,
        name='accuracy',
        marker_color='#EDAE49'
    ))
    fig.add_trace(go.Bar(
        y=df['p1'].values.tolist(),
        x=x,
        name='precision',
        marker_color='#D1495B'
    ))
    fig.add_trace(go.Bar(
        y=df['r1'].values.tolist(),
        x=x,
        name='recall',
        marker_color='#00798C'
    ))
    fig.add_trace(go.Bar(
        y=df['f1'].values.tolist(),
        x=x,
        name='F-score',
        marker_color='#003D5B'
    ))
    fig.update_layout( title='RepeatsDB-lite1', width=1000)
    fig.write_image(cfg.data['plots'] + "/a_p_r_m_dataset_rdb1.png")

    # RDBL 2
    fig2 = go.Figure()

    fig2.add_trace(go.Bar(
        y=df['a2'].values.tolist(),
        x=x,
        name='accuracy',
        marker_color='#EDAE49'
    ))
    fig2.add_trace(go.Bar(
        y=df['p2'].values.tolist(),
        x=x,
        name='precision',
        marker_color='#D1495B'
    ))
    fig2.add_trace(go.Bar(
        y=df['r2'].values.tolist(),
        x=x,
        name='recall',
        marker_color='#00798C'
    ))
    fig2.add_trace(go.Bar(
        y=df['f2'].values.tolist(),
        x=x,
        name='F-score',
        marker_color='#003D5B'
    ))
    fig2.update_layout(title='RepeatsDB-lite2', width=1000)

    fig2.write_image(cfg.data['plots'] + "/a_p_r_m_dataset_rdb2.png")

def plot_metrics_pdbs():
    df = pd.read_csv(cfg.data['data'] + '/analysis_pdbs.csv', sep=',', dtype={'topologies': str})

    # accuracy
    df_plot1 = df[['a1', 'topologies']].rename(columns={"a1": "accuracy"})
    df_plot2 = df[['a2', 'topologies']].rename(columns={"a2": "accuracy"})
    df_plot1['tool'] = 'RepeatsDB-lite1'
    df_plot2['tool'] = 'RepeatsDB-lite2'
    df_final = pd.concat([df_plot1, df_plot2])
    df_final = df_final.sort_values(by=['topologies'], ascending=True)
    df_final[['class','suclass']] = df_final['topologies'].str.split('.',expand=True)


    fig = px.box(df_final, x="topologies", y="accuracy", color='tool')
    fig.update_layout(width=1000)
    fig.write_image(cfg.data['plots'] + "/accuracy_pdbs.png")
    fig = px.box(df_final, x="class", y="accuracy", color='tool')
    fig.update_layout(width=1000)
    fig.write_image(cfg.data['plots'] + "/accuracy_pdbs_classes.png")


    # precision
    df_plot1 = df[['p1', 'topologies']].rename(columns={"p1": "precision"})
    df_plot2 = df[['p2', 'topologies']].rename(columns={"p2": "precision"})
    df_plot1['tool'] = 'RepeatsDB-lite1'
    df_plot2['tool'] = 'RepeatsDB-lite2'
    df_final = pd.concat([df_plot1, df_plot2])
    df_final = df_final.sort_values(by=['topologies'], ascending=True)
    df_final[['class','suclass']] = df_final['topologies'].str.split('.',expand=True)


    fig = px.box(df_final, x="topologies", y="precision", color='tool')
    fig.update_layout(width=1000)
    fig.write_image(cfg.data['plots'] + "/precision_pdbs.png")
    fig = px.box(df_final, x="class", y="precision", color='tool')
    fig.update_layout(width=1000)
    fig.write_image(cfg.data['plots'] + "/precision_pdbs_classes.png")

    # recall
    df_plot1 = df[['r1', 'topologies']].rename(columns={"r1": "recall"})
    df_plot2 = df[['r2', 'topologies']].rename(columns={"r2": "recall"})
    df_plot1['tool'] = 'RepeatsDB-lite1'
    df_plot2['tool'] = 'RepeatsDB-lite2'
    df_final = pd.concat([df_plot1, df_plot2])
    df_final = df_final.sort_values(by=['topologies'], ascending=True)
    df_final[['class','suclass']] = df_final['topologies'].str.split('.',expand=True)


    fig = px.box(df_final, x="topologies", y="recall", color='tool')
    fig.update_layout(width=1000)
    fig.write_image(cfg.data['plots'] + "/recall_pdbs.png")
    fig = px.box(df_final, x="class", y="recall", color='tool')
    fig.update_layout(width=1000)
    fig.write_image(cfg.data['plots'] + "/recall_pdbs_classes.png")


    # recall
    df_plot1 = df[['r1', 'topologies']].rename(columns={"r1": "recall"})
    df_plot2 = df[['r2', 'topologies']].rename(columns={"r2": "recall"})
    df_plot1['tool'] = 'RepeatsDB-lite1'
    df_plot2['tool'] = 'RepeatsDB-lite2'
    df_final = pd.concat([df_plot1, df_plot2])
    df_final = df_final.sort_values(by=['topologies'], ascending=True)
    df_final[['class','suclass']] = df_final['topologies'].str.split('.',expand=True)


    fig = px.box(df_final, x="topologies", y="recall", color='tool')
    fig.update_layout(width=1000)
    fig.write_image(cfg.data['plots'] + "/recall_pdbs.png")
    fig = px.box(df_final, x="class", y="recall", color='tool')
    fig.update_layout(width=1000)
    fig.write_image(cfg.data['plots'] + "/recall_pdbs_classes.png")

    # F score
    df_f1 = df[['f1', 'topologies']].rename(columns={"f1": "f-score"})
    df_f2 = df[['f2', 'topologies']].rename(columns={"f2": "f-score"})
    df_f1['tool'] = 'RepeatsDB-lite1'
    df_f2['tool'] = 'RepeatsDB-lite2'
    df_final = pd.concat([df_f1, df_f2])
    df_final = df_final.sort_values(by=['topologies'], ascending=True)
    df_final[['class','suclass']] = df_final['topologies'].str.split('.',expand=True)

    fig = px.box(df_final, x="topologies", y="f-score", color='tool')
    fig.update_layout(width=1000)
    fig.write_image(cfg.data['plots'] + "/f-score_pdbs.png")
    fig = px.box(df_final, x="class", y="f-score", color='tool')
    fig.update_layout(width=1000)
    fig.write_image(cfg.data['plots'] + "/f-score_pdbs_classes.png")

    # f1 by method
    fig = px.box(df_final, x="tool", y="f-score")
    fig.write_image(cfg.data['plots'] + "/f-score_by_method.png")




# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    # DATASETS
    # create a list of curated pdbs
    # create_input_txt()
    # create_input_negatives()

    # get info such as number of pdb, units...
    # table_class_topology()

    # parsing curated pdbs and predictions into dictionaries have the same format
    # dict_curated = parse_curated_to_dict()
    # dict_predicted1 = parse_predicted1_to_dict()
    # dict_predicted2 = parse_predicted2_to_dict()
    # dict_regions = parse_regions_to_dict()
    # binary file to compare each residue
    dict_to_binary(dict_curated, dict_predicted1, dict_predicted2, dict_regions)
    # add class info
    # classes_dict = curated_to_dict_classes('class')
    # classes_dict2 = curated_to_dict_classes('class_topology')
    # add_class_to_binary(classes_dict, classes_dict2)
    residues_count()



    # STATISTICS
    # # Dataset level
    # calculate balanced accuracy, precision, recall, F score - all residues
    # Barplot
    # calc_fractions_dataset()
    # plot_metrics_dataset()


    # # calculate balanced accuracy, precision, recall, F score - each PDBs
    # calc_fractions_pdbs()
    # plot_metrics_pdbs()
    # plot_f_score()


    # extract_particular_cases()


    # a general  overview
    # calc_fractions()
    # plot_table()
    # matrix plots
    # plot_matrix()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
