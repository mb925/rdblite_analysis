import statistics as stat
import main as functions
import config as cfg
import pandas as pd
import plotly.express as px
from sklearn import preprocessing

def unit_length():
    c = functions.parse_curated_to_dict()
    p1 = functions.parse_predicted1_to_dict()
    p2 = functions.parse_predicted2_to_dict()

    binary = pd.read_csv(cfg.data['data'] + '/binary_pdb.csv', sep=',', dtype={'pdb_residue_id': str})
    data = {}
    binary.groupby('PDB').apply(lambda x: get_intervals(data, c, p1, p2, x))

    df = pd.DataFrame(data).T.reset_index()
    df.columns = ['PDB', 'delta median rdb1', 'delta median rdb2', 'delta std rdb1', 'delta std rdb2']

    # delta median by method

    df_plot1_m = df.rename(columns={"delta median rdb1": "delta median"})
    df_plot2_m = df.rename(columns={"delta median rdb2": "delta median"})

    print(df_plot1_m['delta median'])

    df_plot1_m['tool'] = 'RepeatsDB-lite1'
    df_plot2_m['tool'] = 'RepeatsDB-lite2'
    df_final_m = pd.concat([df_plot1_m, df_plot2_m])



    fig = px.box(df_final_m, x="tool", y="delta median")
    fig.write_image(cfg.data['units'] + "/delta_median_by_method.png")

    # # delta std by method

    df_plot1_s = df.rename(columns={"delta std rdb1": "delta std"})
    df_plot2_s = df.rename(columns={"delta std rdb2": "delta std"})

    df_plot1_s['tool'] = 'RepeatsDB-lite1'
    df_plot2_s['tool'] = 'RepeatsDB-lite2'
    df_final_s = pd.concat([df_plot1_s, df_plot2_s])


    fig = px.box(df_final_s, x="tool", y="delta std")
    fig.write_image(cfg.data['units'] + "/delta_std_by_method.png")

    # print(df)

def get_intervals(data, c, p1, p2, x):
    print(x.iloc[0]['PDB'])
    deltas_c = []
    deltas_p1 = []
    deltas_p2 = []

    intervals_c = c[x.iloc[0]['PDB']]
    for el in intervals_c:

        st_seqres = x.loc[x['pdb_residue_id'] == el[0]]['seqres_index']
        end_seqres = x.loc[x['pdb_residue_id'] == el[1]]['seqres_index']
        if not st_seqres.empty:
            if not end_seqres.empty:
                deltas_c.append(end_seqres.iloc[0] - st_seqres.iloc[0])

    if x.iloc[0]['PDB'] in p1:
        intervals_p1 = p1[x.iloc[0]['PDB']]
        for el in intervals_p1:
            st_seqres = x.loc[x['pdb_residue_id'] == el[0]]['seqres_index']
            end_seqres = x.loc[x['pdb_residue_id'] == el[1]]['seqres_index']
            if not st_seqres.empty:
                if not end_seqres.empty:
                    deltas_p1.append(end_seqres.iloc[0] - st_seqres.iloc[0])

    if x.iloc[0]['PDB'] in p2:
        intervals_p2 = p2[x.iloc[0]['PDB']]
        for el in intervals_p2:
            st_seqres = x.loc[x['pdb_residue_id'] == str(el[0])]['seqres_index']
            end_seqres = x.loc[x['pdb_residue_id'] == str(el[1])]['seqres_index']
            if not st_seqres.empty:
                if not end_seqres.empty:
                    deltas_p2.append(end_seqres.iloc[0] - st_seqres.iloc[0])


    median_c = 0
    median_p1 = 0
    median_p2 = 0
    if len(deltas_c) > 0:
        median_c = round(stat.median(deltas_c), 2)
    if len(deltas_p1) > 0:
        median_p1 = round(stat.median(deltas_p1), 2)
    if len(deltas_p2) > 0:
        median_p2 = round(stat.median(deltas_p2), 2)


    delta_md_rdb1 = abs(median_c - median_p1)
    delta_md_rdb2 = abs(median_c - median_p2)

    std_c = 0
    std_p1 = 0
    std_p2 = 0
    if len(deltas_c) > 1:
        std_c = round(stat.stdev(deltas_c), 2)
    if len(deltas_p1) > 1:
        std_p1 = round(stat.stdev(deltas_p1), 2)
    if len(deltas_p2) > 1:
        std_p2 = round(stat.stdev(deltas_p2), 2)


    delta_std_rdb1 = abs(std_c - std_p1)
    delta_std_rdb2 = abs(std_c - std_p2)

    data[x.iloc[0]['PDB']] = [delta_md_rdb1, delta_md_rdb2, delta_std_rdb1, delta_std_rdb2]


if __name__ == '__main__':

    unit_length()
