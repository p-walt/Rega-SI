import json
import pandas as pd
import csv

file = open(
    '/mnt/c/Users/pcypw1/Documents/PhD/GitHub/ML-for-CH/app/Calculation_data/00c7e83f-2a8a-49ef-aaea-5b64657f912f_results.json')
dic = json.load(file)

with open('/mnt/c/Users/pcypw1/Documents/PhD/GitHub/ML-for-CH/app/Calculation_data/00c7e83f-2a8a-49ef-aaea-5b64657f912f_results.csv', 'w') as file:
    header = ['compound_name', 'smiles', 'radical', 'Site', 'hf_reactant_energy', 'hf_ts_energy','hf_activation_energy', 'hf_ratio', 'Fukui Fa+', 'Fukui Fa-', 'Fukui Fa0']
    writer = csv.writer(file)
    writer.writerow(header)
    for com in dic.keys():
        for radical in dic[com]['radicals'].values():
            for key, value in radical.items():
                if str.isdigit(key):
                    row = [dic[com]['compound_name'], dic[com]['smiles'], radical['radical'], key, radical['hf_reactant_energy'], radical[key]['hf_ts_energy'], radical[key]['hf_act_energy'], radical[key]['hf_ratio'], radical[key]['fa+'], radical[key]['fa-'], radical[key]['fa0']]
                    writer.writerow(row)
            # print(radical['hf_reactant_energy'], dic['Name'], radical['Dept num'], radical['Dept name'], item['First name'], item['Last name'])