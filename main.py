# This is a sample Python script.
import ast

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.
import pandas as pd
import json
import config as cfg

# create input.txt file from entries
def create_input_txt():
    f = open(cfg.data['data'] + '/entries.json')
    data = json.load(f)
    pdb_list = []
    for i in data:

        if i['reviewed'] == True:
            print(i['reviewed'])
            pdb_list.append(i['repeatsdb_id'])

    print(pdb_list)
    pdb_list = set(pdb_list)

    with open(cfg.data['data'] + '/input.txt', 'a') as the_file:
        for el in pdb_list:
            the_file.write(el[0:4] + ' ' + el[4] + '\n')




# calculates intersections among 2 intervals
def intersections(a,b):
    ranges = []
    i = j = 0
    while i < len(a) and j < len(b):
        a_left, a_right = a[i]
        b_left, b_right = b[j]

        if a_right < b_right:
            i += 1
        else:
            j += 1

        if a_right >= b_left and b_right >= a_left:
            end_pts = sorted([a_left, a_right, b_left, b_right])
            middle = [end_pts[1], end_pts[2]]
            ranges.append(middle)

    ri = 0
    while ri < len(ranges)-1:
        if ranges[ri][1] == ranges[ri+1][0]:
            ranges[ri:ri+2] = [[ranges[ri][0], ranges[ri+1][1]]]

        ri += 1
    # from ranges need to calculate and sum the deltas
    return ranges


def parse_curated_to_dict():
    print('entries dict')
    dict_curated = {}
    f = open(cfg.data['data'] + 'entries.json')
    data = json.load(f)
    for i in data:
        if i["repeatsdb_id"] in dict_curated:
            dict_curated[i["repeatsdb_id"]].append([int(i['start']), int(i['end'])])
        if i["repeatsdb_id"] not in dict_curated:
            dict_curated[i["repeatsdb_id"]] = []
            dict_curated[i["repeatsdb_id"]].append([int(i['start']), int(i['end'])])
    return dict_curated

def parse_predicted_to_dict():
    print('predicted dict')
    dict_predicted = {}
    data = pd.read_csv(cfg.data['data'] + '/results.csv', sep=',')

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

def calculate_overlap(curated, predicted):
    with open(cfg.data['data'] + '/analysis.csv', 'a') as the_file:
        the_file.write('PDB' + ',' + 'OVERLAP' + '\n')
        for key in curated:
            if key in predicted:
                overlap = intersections(curated[key], predicted[key])
                print(overlap)
                residues_overlap = count_residues(overlap)
                residues_curated = count_residues(curated[key])
                percentage = 2 * residues_overlap/(residues_curated + residues_overlap)
                the_file.write(key + ',' + str(round(percentage, 2)) + '\n')




def count_residues(deltas):
    d = 0
    for delta in deltas:
        d += delta[1] - delta[0]
    return d

def average():
    data = pd.read_csv(cfg.data['data'] + '/analysis.csv', sep=',')
    mean = data["OVERLAP"].mean()
    print(mean)

def curated_to_dict_classes(level):
    print('entries dict')
    dict_curated = {}
    f = open(cfg.data['data'] + '/entries.json')
    data = json.load(f)
    for i in data:
        if i["repeatsdb_id"] in dict_curated:
            dict_curated[i["repeatsdb_id"]].append(i[level])
        if i["repeatsdb_id"] not in dict_curated:
            dict_curated[i["repeatsdb_id"]] = []
            dict_curated[i["repeatsdb_id"]].append(i[level])
    return dict_curated

def predicted_to_dict_classes(level):
    print('predicted dict')
    dict_predicted = {}
    data = pd.read_csv(cfg.data['data'] + '/results.csv', sep=',')

    ontology = pd.read_csv(cfg.data['data'] + '/ontology.csv', sep=',')
    for i in data.iterrows():
        if i[1][2] == 'no regions':
            continue
        cl = i[1][level] # convert in number
        print(cl)
        cl_name = ontology.loc[ontology['Name'].str.contains(cl)]['Code'].to_string(index=False)
        print(cl_name)


        if i[1][0] in dict_predicted:
            dict_predicted[i[1][0]].append(cl_name)
        if i[1][0] not in dict_predicted:
            dict_predicted[i[1][0]] = []
            dict_predicted[i[1][0]].append(cl_name)
    return dict_predicted


def compare_classes(curated, predicted):
    with open(cfg.data['data'] + '/analysis_classes.csv', 'a') as the_file:
        for key in curated:
            if key in predicted:
                # print(len(set(curated[key]) & set(predicted[key])))
                the_file.write(key + ',' + str(len(set(curated[key]) & set(predicted[key]))) + '\n')
                # print(list(set(set(curated[key])).intersection(set(predicted[key]))))

def count_classes(level):
    data = pd.read_csv(cfg.data['data'] + '/analysis_classes.csv', sep=',')
    # print(data.loc[data["CLASS"] == 1])
    count = len(data.loc[data[level] >= 1]) / len(data)
    print(level, count)

def classes_analysis():
    with open(cfg.data['data'] + '/analysis_classes.csv', 'w') as the_file:
        # the_file.write('PDB' + ',' + 'CLASS' + ',' + 'TOPOLOGY' + ',' + 'FOLD' + '\n')
        the_file.write('PDB' + ',' + 'TOPOLOGY' + '\n')
        # check class
    # dict_curated_classes = curated_to_dict_classes('class')
    # dict_predicted_classes = predicted_to_dict_classes(4)
    # compare_classes(dict_curated_classes, dict_predicted_classes)
    # count_classes('CLASS')
    a = curated_to_dict_classes('class_topology')
    b = predicted_to_dict_classes(5)
    compare_classes(a, b)
    count_classes('TOPOLOGY')


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    create_input_txt()

    # check overlap
    # dict_curated = parse_curated_to_dict()
    # dict_predicted = parse_predicted_to_dict()
    #
    # calculate_overlap(dict_curated, dict_predicted)
    # average()

    # classes_analysis()



# See PyCharm help at https://www.jetbrains.com/help/pycharm/
