import uuid
import os
import ast
from calculations import calculations_bp
from flask_login import login_required
from flask import render_template, session
from app import db, pycharm
from flask_login import current_user
from pony.orm import select
from operator import itemgetter
from func_calc import am1_calcs
from auxiliary import get_image_data
from pathlib import Path


@calculations_bp.route('/my_calculations', methods=['GET', 'POST'])
@login_required
def my_calculations():
    # TODO RA: Docstring
    calcs = []
    for calc in select(u for u in db.Calculation if current_user.email == u.email):
        calc = calc.to_dict()
        if len(calc["smiles"]) == 1:
            calc["mode"] = "single"
        else:
            calc["mode"] = "multiple"
        calc["smiles"] = ", ".join(calc["smiles"])
        if len(calc["smiles"]) > 20:
            calc["smiles"] = calc["smiles"][0:19] + "..."
        calc["radical"] = ", ".join(calc["radical"])
        if len(calc["radical"]) > 20:
            calc["radical"] = calc["radical"][0:16] + "..."
        calc["name"] = ", ".join(calc["name"])
        if len(calc["name"]) > 20:
            calc["name"] = calc["name"][0:16] + "..."
        calc["date_time"] = calc["date_time"][:-7]
        calcs.append(calc)
    calcs = sorted(calcs, key=itemgetter('date_time'))
    return render_template("my_calculations.html", calcs=calcs)


@calculations_bp.route('/reload/<data>', methods=['GET', 'POST'])
@login_required
def reload(data):
    # TODO RA: Docstring

    data = ast.literal_eval(data)
    if pycharm == 1:
        if "Calculation_data" not in os.getcwd():
            start = Path.cwd() / 'ML-for-CH' / 'app' / 'Calculation_data'
    else:
        if "Calculation_data" not in os.getcwd():
            start = Path.cwd() / 'Calculation_data'
    ran_string = str(uuid.uuid4())
    session["random_name"] = ran_string
    molecule_dir = str(start / ran_string)
    # if single
    if data["mode"] == "single":
        results_dict, molecule_sites = am1_calcs.calculate_am1_energies_single(data["name"], data["smiles"],
                                                                               data["radical"], start, molecule_dir)
        compounds = results_dict.keys()
        row_headers = list(list(results_dict.values())[0].values())[2].keys()
        compound_eas = {}
        compound_sites = {}
        compound_ratios = {}
        radicals = {}
        for compound in compounds:
            activation_energies = {}
            ratios = {}
            radicals[compound] = data["radical"]
            for site in molecule_sites:
                activation_energies[str(site)] = (
                results_dict[compound]['radicals'][str(data["radical"])][str(site)]["am1_act_energy"])
                ratios[str(site)] = (results_dict[compound]['radicals'][str(data["radical"])][str(site)]["am1_ratio"])
            compound_eas[compound] = activation_energies
            compound_ratios[compound] = ratios
            compound_sites[compound] = molecule_sites
        image_data = get_image_data(data["smiles"], molecule_sites)
        return render_template("complete.html", image_data=[image_data], compounds=compounds, row_headers=row_headers,
                               results_dict=results_dict, compound_eas=compound_eas, compound_sites=compound_sites,
                               compound_ratios=compound_ratios, random_name=session["random_name"], mode="single",
                               radicals=radicals)
    # else multiple
    else:
        new_dir = str(start / data["random_string"])
        results_dict, sites_dict, smiles_dict, radicals_dict = am1_calcs.calculate_am1_energies_list(
            Path(f'{new_dir}.csv'), start, Path(f'{new_dir}'), ran_string)
        compounds = results_dict.keys()
        row_headers = list(list(results_dict.values())[0].values())[2].keys()
        smiles_list = smiles_dict.keys()
        compound_eas = {}
        compound_sites = {}
        compound_ratios = {}
        radicals_list = []
        for smile in smiles_list:
            radicals_list.append(radicals_dict[smile])
            activation_energies = {}
            ratios = {}
            for site in sites_dict[smile]:
                activation_energies[str(site)] = (
                results_dict[smiles_dict[smile]]['radicals'][radicals_dict[smile]][str(site)]["am1_act_energy"])
                ratios[str(site)] = (
                results_dict[smiles_dict[smile]]['radicals'][radicals_dict[smile]][str(site)]["am1_ratio"])
            compound_eas[smiles_dict[smile]] = activation_energies
            compound_ratios[smiles_dict[smile]] = ratios
            compound_sites[smiles_dict[smile]] = sites_dict[smile]
        radicals = {}
        for i in range(len(compounds)):
            radicals[list(compounds)[i]] = radicals_list[i]
        image_data = []
        smiles_list = list(smiles_list)
        for i in range(len(smiles_list)):
            image_data.append(get_image_data(smiles_list[i], sites_dict[smiles_list[i]]))
        return render_template("complete.html", image_data=image_data, compounds=compounds, row_headers=row_headers,
                               results_dict=results_dict, compound_eas=compound_eas, compound_sites=compound_sites,
                               compound_ratios=compound_ratios, random_name=session["random_name"], mode="multiple",
                               radicals=radicals)
