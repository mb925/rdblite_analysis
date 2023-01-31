import main as m
import json
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import config as cfg
from Bio import SeqIO
from pymongo import MongoClient

def count_not_trp():

    rdb2 = pd.read_csv(cfg.data['data'] + '/rpdblite2_predictions.csv', sep=',')
    topologies = pd.read_csv(cfg.data['data'] + '/analysis_pdbs.csv', sep=',', dtype={'topologies': str})[['PDB', 'topologies']]

    f1 = pd.read_csv(cfg.data['data'] + '/particular_cases/rdb1_f1_0.csv', sep=',')
    f1_merged = rdb2.merge(f1, how='right', on='PDB')
    f1_merged_all = topologies.merge(rdb2, how='left', on='PDB')

    # rdblite2 - F score 0
    df2 = pd.DataFrame()
    detected2 = []
    for i in f1_merged.iterrows():
        if i[1][2] == 'no regions':
            detected2.append('Not detected')
        else:
            detected2.append('Detected')
    df2['Detection'] = detected2

    fig1 = go.Figure(data =[go.Table(header=dict(values=['Tool','Detected', 'Not detected', 'Total PDBs']),
    cells = dict(values=[['Repeatsdb-lite2'],[len(df2.loc[df2['Detection'] == 'Detected'])],
                         [len(df2.loc[df2['Detection'] == 'Not detected'])],[len(df2)]]

                 ))
    ])

    fig1.update_layout(title="RDB2 F-score = 0")
    fig1.write_image(cfg.data['plots'] + "/particular_cases/rdb2_f1_0.png")

    # rdblite2 - F score all

    detected2_all = []
    for i in f1_merged_all.iterrows():
        if i[1][2] == 'no regions':
            detected2_all.append('Not detected')
        else:
            detected2_all.append('Detected')
    f1_merged_all['Detection'] = detected2_all

    topologies_sorted = f1_merged_all.sort_values(by=['topologies'], ascending=True)
    counts = topologies_sorted.groupby(['topologies', 'Detection'])['topologies'].size().reset_index(name='counts')


    fig1 = px.bar(counts, x="topologies", y="counts", color="Detection", title="Detection level by topology - Repeatsdb-lite2", barmode='stack')
    fig1.write_image(cfg.data['plots'] + "/particular_cases/detection-topologies-Repeatsdb-lite2.png")





def classes_analysis():
    with open(cfg.data['data'] + '/analysis_classes.csv', 'w') as the_file:
        the_file.write('PDB' + ',' + 'CLASS' + ',' + 'TOPOLOGY' + ',' + 'FOLD' + '\n')

    a = curated_to_dict_classes('class')
    b = predicted_to_dict_classes('class')
    compare_classes(a, b)


def plot_classes():
    rdb2 = pd.read_csv(cfg.data['data'] + '/particular_cases/classes-rdb2.csv', sep=',')
    fig1 = px.histogram(rdb2, x=rdb2["CLASS"], title="classes - Repeatsdb-lite2")
    fig1.write_image(cfg.data['plots'] + "/particular_cases/classes - Repeatsdb-lite2.png")


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


def predicted_to_dict_classes(level):
    print('predicted dict')
    dict_predicted = {}
    rdb2 = pd.read_csv(cfg.data['data'] + '/rpdblite2_predictions.csv', sep=',')

    mcc = pd.read_csv(cfg.data['data'] + '/particular_cases/mcc_0.csv', sep=',')
    mcc_merged = rdb2.merge(mcc, how='right', on='PDB')

    ontology = pd.read_csv(cfg.data['data'] + '/ontology.csv', sep=',')
    for i in mcc_merged.iterrows():
        if i[1][2] == 'no regions':
            continue
        cl = i[1][level]
        print('cl', cl)

        # convert in number
        cl_num = ontology.loc[ontology['Name'].str.contains(cl)]['Code'].to_string(index=False)
        print('num', cl_num)

        if i[1][0] in dict_predicted:
            dict_predicted[i[1][0]].append(cl_num)
        if i[1][0] not in dict_predicted:
            dict_predicted[i[1][0]] = []
            dict_predicted[i[1][0]].append(cl_num)
    return dict_predicted

def compare_classes(curated, predicted):
    with open(cfg.data['data'] + '/particular_cases/classes-rdb2.csv', 'a') as the_file:
        the_file.write('PDB' + ',' + 'CLASS' + '\n')

        for key in curated:
            if key in predicted:
                # print(len(set(curated[key]) & set(predicted[key])))
                the_file.write(key + ',' + str(len(set(curated[key]) & set(predicted[key]))) + '\n')
                # print(list(set(set(curated[key])).intersection(set(predicted[key]))))


def dict_to_binary(d_pred1, d_pred2):
    # EBI file parsing
    df = pd.DataFrame()

    list_negatives = pd.read_csv(cfg.data['data'] + '/input_negatives.txt', sep=' ', names=["PDB", "CHAIN"])

    # Connect to mongo (BioDBs seqres)
    client = MongoClient(cfg.db_host, cfg.db_port)
    db = client.biodbs

    for pdb in list_negatives.iterrows():
        l = m.find_seqres(db, pdb[1]['PDB'], pdb[1]['CHAIN'])
        pdb = pdb[1]['PDB'] + pdb[1]['CHAIN']
        df_residues = pd.DataFrame.from_records(l)
        if df_residues.empty:
            print('empty ', pdb)
            continue
        residues = df_residues.sort_values(by=['seqres_index'], ascending=True)[
            ['pdb_id', 'pdb_chain', 'seqres_index', 'pdb_residue_id']]

        residues['CURATED'] = 0
        residues['RDB1'] = 0
        residues['RDB2'] = 0

        if pdb in d_pred1:
            for interval1 in d_pred1[pdb]:
                m.fill_repeated_interval(pdb, residues, interval1, 'RDB1')
        if pdb in d_pred2:
            for interval2 in d_pred2[pdb]:
                m.fill_repeated_interval(pdb, residues, interval2, 'RDB2')

        df = pd.concat([df, residues])


    df["PDB"] = df['pdb_id'] + df["pdb_chain"]
    df.to_csv(cfg.data['cases'] + '/binary_pdb_negatives.csv', index=False)


def create_table():
    data = []
    df = pd.read_csv(cfg.data['cases'] + '/binary_pdb_negatives.csv')
    df_tn1 = len(df.loc[df['RDB1'] == 0])
    df_fp1 = len(df.loc[df['RDB1'] == 1])
    total_1 = df_tn1 + df_fp1
    df_tn2 = len(df.loc[df['RDB2'] == 0])
    df_fp2 = len(df.loc[df['RDB2'] == 1])
    total_2 = df_tn2 + df_fp2
    data.append(['RDBLITE1', df_tn1, df_fp1, total_1])
    data.append(['RDBLITE2', df_tn2, df_fp2, total_2])
    table = pd.DataFrame(data, columns=['TOOL', 'TN', 'TP', 'TOTAL'])
    table.to_csv(cfg.data['cases'] + '/table_negatives.csv', index=False)
    print(df_fp1)

if __name__ == '__main__':

    # count_not_trp()
    # classes_analysis()
    # plot_classes()

    ## Negatives
    # data = pd.read_csv(cfg.data['data'] + '/particular_cases/input_negatives.csv', sep=' ')
    # data = data['PDB']
    # data.to_csv(cfg.data['data'] + '/particular_cases/list_Ebi.csv', index=False)
    # parse_fasta(cfg.data['data'] + '/particular_cases/ebi.fasta')
    # dict_predicted1 = m.parse_predicted1_to_dict(cfg.data['data'] + '/rpdblite1_predictions.tsv')
    # dict_predicted2 = m.parse_predicted2_to_dict(cfg.data['data'] + '/rpdblite2_predictions.csv')
    # dict_to_binary(dict_predicted1, dict_predicted2)
    create_table()
