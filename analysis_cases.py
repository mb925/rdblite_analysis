import sys
import json
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import config as cfg

def count_not_trp():

    rdb2 = pd.read_csv(cfg.data['data'] + '/rpdblite2_predictions.csv', sep=',')
    rdb1 = pd.read_csv(cfg.data['data'] + '/rpdblite1_predictions.tsv', sep='\t', names=["PDB", "TYPE", "START", "END"])
    topologies = pd.read_csv(cfg.data['data'] + '/analysis_pdbs.csv', sep=',', dtype={'topologies': str})[['PDB', 'topologies']]

    mcc = pd.read_csv(cfg.data['data'] + '/particular_cases/mcc_0.csv', sep=',')
    mcc_merged = rdb2.merge(mcc, how='right', on='PDB')

    # rdblite2
    df2 = pd.DataFrame()
    detected2 = []
    for i in mcc_merged.iterrows():
        if i[1][2] == 'no regions':
            detected2.append('Not detected')
        else:
            detected2.append('Detected')
    df2['Detection'] = detected2

    fig2 = go.Figure()
    fig2.add_trace(go.Bar(
        y=['Tandem repeat'],
        x=[len(df2.loc[df2['Detection'] == 'Detected'])],
        name='Detected',
        orientation='h',
        marker=dict(
            color='rgba(246, 78, 139, 0.6)',
            line=dict(color='rgba(246, 78, 139, 1.0)', width=3)
        )
    ))
    fig2.add_trace(go.Bar(
        y=['Tandem repeat'],
        x=[len(df2.loc[df2['Detection'] == 'Not detected'])],
        name='Not detected',
        orientation='h',
        marker=dict(
            color='rgba(58, 71, 80, 0.6)',
            line=dict(color='rgba(58, 71, 80, 1.0)', width=3)
        )
    ))

    fig2.update_layout(barmode='stack', height=300, title="mcc 0 - Repeatsdb-lite2")
    fig2.write_image(cfg.data['plots'] + "/particular_cases/mcc_0 - Repeatsdb-lite2.png")

    # rdblite1
    df1 = pd.DataFrame()

    detected1 = []

    rdb1_values = []
    rdb1_pdbs = []
    for i in rdb1.iterrows():
        pdb = i[1]['PDB'].split("_")[0] + i[1]['PDB'].split("_")[1]
        pdb = pdb[3:]
        rdb1_values.append(pdb)

    for mc in mcc.iterrows():
        if mc[1]['PDB'] in rdb1_values:
            rdb1_pdbs.append(mc[1]['PDB'])
            detected1.append('Detected')
        else:
            rdb1_pdbs.append(mc[1]['PDB'])
            detected1.append('Not detected')

    df1['Detection'] = detected1
    df1['PDB'] = rdb1_pdbs
    df1_classes = topologies.merge(df1, how='right', on='PDB')

    topologies_sorted = df1_classes.sort_values(by=['topologies'], ascending=True)
    topologies = list(topologies_sorted['topologies'].drop_duplicates().values)

    counts = topologies_sorted.groupby(['topologies', 'Detection'])['topologies'].size().reset_index(name='counts')

    fig1 = go.Figure(data=[
        go.Bar(name='Not detected', x=topologies, y=list(counts.loc[counts['Detection'] == 'Not detected']['counts'])),
        go.Bar(name='Detected', x=topologies, y=list(counts.loc[counts['Detection'] == 'Detected']['counts']))
    ])
    # Change the bar mode
    fig1.update_layout(barmode='stack', title="mcc 0 - Repeatsdb-lite1")
    fig1.write_image(cfg.data['plots'] + "/particular_cases/mcc_0 - topologies - Repeatsdb-lite1.png")

    fig1 = go.Figure()
    fig1.add_trace(go.Bar(
        y=['Tandem repeat units'],
        x=[len(df1.loc[df1['Detection'] == 'Detected'])],
        name='Detected',
        orientation='h',
        marker=dict(
            color='rgba(246, 78, 139, 0.6)',
            line=dict(color='rgba(246, 78, 139, 1.0)', width=3)
        )
    ))
    fig1.add_trace(go.Bar(
        y=['Tandem repeat units'],
        x=[len(df1.loc[df1['Detection'] == 'Not detected'])],
        name='Not detected',
        orientation='h',
        marker=dict(
            color='rgba(58, 71, 80, 0.6)',
            line=dict(color='rgba(58, 71, 80, 1.0)', width=3)
        )
    ))

    fig1.update_layout(barmode='stack', height=300, title="mcc 0 - Repeatsdb-lite1")
    fig1.write_image(cfg.data['plots'] + "/particular_cases/mcc_0 - Repeatsdb-lite1.png")





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


def mcc_by_class():
    tools = ['RepeatsDB-lite 1', 'RepeatsDB-lite 2']

    fig = go.Figure(data=[
        go.Bar(name='SF Zoo', x=tools, y=[20, 14]),
        go.Bar(name='LA Zoo', x=tools, y=[12, 18])
    ])
    # Change the bar mode
    fig.update_layout(barmode='stack')
    fig.show()

if __name__ == '__main__':

    count_not_trp()
    # mcc_by_class()
    # classes_analysis()
    # plot_classes()
