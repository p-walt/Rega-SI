import os
from app import pycharm, db
from main import main_bp
from flask import render_template, request, url_for, redirect, send_file, session
from flask_login import current_user
import pandas as pd
from io import BytesIO
from func_calc import am1_calcs
from auxiliary import valid_smiles, valid_smiles_multiple, check_size, check_size_multiple, get_image_data
from forms import SubmitForm
import uuid
from pathlib import Path
from datetime import datetime


# Render main page
@main_bp.route('/home', methods=['GET', 'POST'])
@main_bp.route('/', methods=['GET', 'POST'])
def home():
    """
    Handles the home page, accessible via the root endpoints. This view function
    supports both GET and POST HTTP methods. On a POST request, it redirects the
    user to appropriate routes based on the selection input. If the user selects
    "Input Molecule", they are redirected to a single molecule input form route.
    If the user selects "Upload File", they are redirected to a multiple molecules
    upload file route. On a GET request, it renders the homepage template.

    :return: A rendered homepage template for GET requests or a redirect response for valid POST requests.
    :rtype: str or flask.Response
    """
    if request.method == "POST":
        if request.form["choice"] == "Input Molecule":
            return redirect(url_for("main.single"))
        if request.form["choice"] == "Upload File":
            return redirect(url_for("main.multiple"))
    return render_template("home.html")


@main_bp.route('/single', methods=['GET', 'POST'])
def single():
    """
    Handles interactions for submitting a single molecule job in a web-based application. Provides UI for
    input fields where users can enter molecule details such as name, SMILES string, and radical type,
    validates the input, calculates molecule properties, and generates results for rendering on success.
    Implements GET and POST HTTP methods. If POST and form submission are successful, redirects to display
    detailed computation results.

    :param pycharm: Indicator to determine the root directory for 'Calculation_data'.
    :type pycharm: int
    :param start: The directory path where calculated data will be stored.
    :type start: Path
    :param errors: List of validation error messages, if any occurred during form submission.
    :type errors: list
    :param page: An instance of the molecule submission form.
    :type page: SubmitForm
    :param ran_string: A unique identifier for storing session data and results.
    :type ran_string: str
    :param molecule_dir: Directory path where computation files for the molecule will be stored.
    :type molecule_dir: str
    :param results_dict: Dictionary containing computation results for submitted molecules.
    :type results_dict: dict
    :param molecule_sites: List of reaction sites on the molecule being calculated.
    :type molecule_sites: list
    :param compounds: Keys representing individual compound results from `results_dict`.
    :type compounds: Iterable
    :param row_headers: A list of header fields for result rows derived from the results_dict.
    :type row_headers: list
    :param compound_eas: Dictionary mapping each compound to its activation energies across sites.
    :type compound_eas: dict
    :param compound_sites: Dictionary mapping each compound to its respective molecule sites.
    :type compound_sites: dict
    :param compound_ratios: Dictionary storing computed molecular reaction ratios across active sites.
    :type compound_ratios: dict
    :param radicals: Mapping of each compound to the radical provided for computations.
    :type radicals: dict
    :param image_data: Encoded string data representing the visual molecule depiction.
    :type image_data: str
    :param dateTimeObj: The current timestamp used for user-specific logging.
    :type dateTimeObj: datetime

    :raises ValidationError: If the molecule data fields (e.g., name, SMILES, radical type) are not provided
                              or contain invalid data.
    :raises SizeLimitError: If the molecule exceeds the maximum allowable size for calculations.
    :raises InvalidSMILES: When the provided SMILES string fails syntax or structural validation checks.

    :return: On successful form submission and computations, renders the result page, else renders the input form
             page with error messages.
    :rtype: werkzeug.wrappers.response.Response
    """
    if pycharm == 1:
        if "Calculation_data" not in os.getcwd():
            start = Path.cwd() / 'Calculation_data'
    else:
        if "Calculation_data" not in os.getcwd():
            start = Path.cwd() / 'Calculation_data'
    errors = []
    page = SubmitForm()
    # submitting a job
    if page.submit_single.data:
        if not page.name.data:
            errors.append("Please choose a name")
            return render_template('single.html', form=page, errors=errors)
        if not page.smiles.data:
            errors.append("Please enter a SMILES string")
            return render_template('single.html', form=page, errors=errors)
        if page.radical.data == "-select-":
            errors.append("Please choose a radical")
            return render_template('single.html', form=page, errors=errors)
        size_check = check_size(page.smiles.data)
        if size_check > 100:
            errors.append("The maximum size for a molecule is 100 non-hydrogen atoms")
            return render_template('single.html', form=page, errors=errors)
        if valid_smiles(page.smiles.data) == 0:
            errors.append(f"The SMILES you entered: {page.smiles.data} is not valid, please try another molecule")
            return render_template('single.html', form=page, errors=errors)

        ran_string = str(uuid.uuid4())
        session["random_name"] = ran_string
        molecule_dir = str(start / ran_string)
        print(f'HPC CALCS ARE {page.hpc_calcs.data}')
        results_dict, molecule_sites = am1_calcs.calculate_am1_energies_single(page.name.data, page.smiles.data,
                                                                               page.radical.data, start, molecule_dir,
                                                                               hpc_calcs=page.hpc_calcs.data)
        compounds = results_dict.keys()
        row_headers = list(list(results_dict.values())[0].values())[2].keys()
        compound_eas = {}
        compound_sites = {}
        compound_ratios = {}
        radicals = {}
        for compound in compounds:
            activation_energies = {}
            ratios = {}
            radicals[compound] = page.radical.data
            if not molecule_sites:
                activation_energies['No Sites of Reaction Found'] = 'No Sites of Reaction Found'
                ratios['No Sites of Reaction Found'] = 'No Sites of Reaction Found'
            for site in molecule_sites:
                if page.hpc_calcs.data:
                    activation_energies[str(site)] = (
                        results_dict[compound]['radicals'][str(page.radical.data)][str(site)]["hf_act_energy"])
                    ratios[str(site)] = (
                        results_dict[compound]['radicals'][str(page.radical.data)][str(site)]["hf_ratio"])
                else:
                    activation_energies[str(site)] = (
                        results_dict[compound]['radicals'][str(page.radical.data)][str(site)]["am1_act_energy"])
                    ratios[str(site)] = (
                        results_dict[compound]['radicals'][str(page.radical.data)][str(site)]["am1_ratio"])
            compound_eas[compound] = activation_energies
            compound_ratios[compound] = ratios
            compound_sites[compound] = molecule_sites
        image_data = get_image_data(page.smiles.data, molecule_sites)
        if current_user.is_authenticated:
            dateTimeObj = datetime.now()
            db.Calculation(email=current_user.email, name=[page.name.data], smiles=[page.smiles.data],
                           radical=[page.radical.data],
                           date_time=str(dateTimeObj), mode="single", random_string=ran_string)
        return render_template("complete.html", image_data=[image_data], compounds=compounds, row_headers=row_headers,
                               results_dict=results_dict, compound_eas=compound_eas, compound_sites=compound_sites,
                               compound_ratios=compound_ratios, random_name=session["random_name"], mode="single",
                               radicals=radicals)
    return render_template("single.html", form=page, errors=errors)


@main_bp.route('/multiple', methods=['GET', 'POST'])
def multiple():
    """
    Handles the logic for the '/multiple' route in the application. This endpoint allows users to upload a
    file containing molecular data and validate, process, and compute results based on the molecules provided.
    The calculations may vary depending on the use of HPC mode or standard AM1 calculations.

    The function checks for the validity of the uploaded file, processes it for further computations, and
    provides the user with results. It is sensitive to the format and size of the uploaded file, ensuring
    strict validity checks and error feedback. Successful submissions trigger AM1 energy calculations for
    molecules in the file. The results from the computations, alongside supplementary visualizations,
    are displayed to the user.

    :error: Renders error messages for various validation failures such as an invalid file type, excessive
            molecule size, or invalid SMILES strings in the input file.
    :param: None
    :returns: HTML template with results, visualized data, or errors for end users.
    :rtype: str
    """
    if pycharm == 1:
        if "Calculation_data" not in os.getcwd():
            start = Path.cwd() / 'ML-for-CH' / 'app' / 'Calculation_data'
    else:
        if "Calculation_data" not in os.getcwd():
            start = Path.cwd() / 'Calculation_data'
    errors = []
    page = SubmitForm()
    print('AT START')
    print(page.submit_multiple.data)
    if page.submit_multiple.data:
        print('SUBMITTED')
        # if not page.hpc_calcs.data:
        #     errors.append("Please tick box")
        #     return render_template('multiple.html', form=page, errors=errors)
        if not page.file.data:
            errors.append("Please choose a file")
            return render_template('multiple.html', form=page, errors=errors)
        try:
            f = BytesIO(page.file.data.read())
            f.seek(0)
            f = pd.read_csv(f)
            size_check = check_size_multiple(f)
            if size_check > 100:
                errors.append("The maximum size for a molecule is 100 non-hydrogen atoms")
                return render_template('multiple.html', form=page, errors=errors)
            valid_smiles_m_result, smiles_culprit = valid_smiles_multiple(f)
            if valid_smiles_m_result is False:
                errors.append(
                    f"One of the SMILES you entered: {smiles_culprit} is not valid, please try another molecule")
                return render_template('multiple.html', form=page, errors=errors)
            ran_string = str(uuid.uuid4())
            session["random_name"] = ran_string
            new_dir = str(start / ran_string)
            f.to_csv(Path(f'{new_dir}.csv'), index=False)
            os.mkdir(Path(f'{new_dir}'))
        except:
            errors.append("Check file is .csv and try again")
            return render_template('multiple.html', form=page, errors=errors)
        print(f'HPC CALCS ARE {page.hpc_calcs.data}')
        results_dict, sites_dict, smiles_dict, radicals_dict = am1_calcs.calculate_am1_energies_list(
            Path(f'{new_dir}.csv'), start, Path(f'{new_dir}'), ran_string, hpc_calcs=page.hpc_calcs.data)
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
            print(sites_dict[smile])
            if not sites_dict[smile]:
                activation_energies['No Sites of Reaction Found'] = 'No Sites of Reaction Found'
                ratios['No Sites of Reaction Found'] = 'No Sites of Reaction Found'
            for site in sites_dict[smile]:
                if page.hpc_calcs.data:
                    activation_energies[str(site)] = (
                        results_dict[smiles_dict[smile]]['radicals'][radicals_dict[smile]][str(site)]["hf_act_energy"])
                    ratios[str(site)] = (
                        results_dict[smiles_dict[smile]]['radicals'][radicals_dict[smile]][str(site)]["hf_ratio"])
                else:
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
        if current_user.is_authenticated:
            dateTimeObj = datetime.now()
            db.Calculation(email=current_user.email, name=list(compounds), smiles=smiles_list, radical=radicals_list,
                           date_time=str(dateTimeObj), mode="multiple", random_string=ran_string)
        return render_template("complete.html", image_data=image_data, compounds=compounds, row_headers=row_headers,
                               results_dict=results_dict, compound_eas=compound_eas, compound_sites=compound_sites,
                               compound_ratios=compound_ratios, random_name=session["random_name"], mode="multiple",
                               radicals=radicals)
    return render_template('multiple.html', form=page, errors=errors)


@main_bp.route('/download_results', methods=['GET', 'POST'])
def download_results():
    """
    Handles the download of calculation result files by serving them as an attachment
    to the client. The file is identified through a random name provided as a query
    parameter, and its path is determined based on the application's directory
    structure and execution context.

    :param random_name: Query parameter received from the client, representing
                        the random name of the file to be downloaded.
    :type random_name: str

    :raises FileNotFoundError: If the specified file does not exist in the
                               constructed file path.
    :raises TypeError: If the provided random_name is of an invalid type.

    :return: A response object that serves the specified file as an attachment to
             the client.
    :rtype: werkzeug.wrappers.Response
    """
    if pycharm == 1:
        if "Calculation_data" not in os.getcwd():
            start = Path.cwd() / 'ML-for-CH' / 'app' / 'Calculation_data'
    else:
        if "Calculation_data" not in os.getcwd():
            start = Path.cwd() / 'Calculation_data'
    random_name = request.args.get('random_name')
    return send_file(Path(f'{start}/{random_name}_results.json'), as_attachment=True)
