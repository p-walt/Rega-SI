 
Rega is an automated high-throughput transition state location program for the generation of machine learning datasets to predict regioselectivity. While the current state of the Rega is designed for the calculation of C-H functionalisation transition states, its modularity allows for easy adaptation to any system of interest.  
The basic functionality of Rega is outlined below:

•	Identification of potential sites of reaction within input molecules

•	Initial semi-empirical transition state search performed on the local machine

•	Ab-initio calculation setup and submission to HPC facility of choice

•	Monitoring of remote HPC calculation and active correction of failing calculations to maximise success rate.

In order to run the graphical user interface (GUI), the Python package Flask is used which parses the user’s inputs and feeds the required information into Rega. It is therefore preferred to run Rega through an IDE such as PyCharm with the run configuration set to a Flask server with the startup script set to app/app.py.
The python environment required to run Rega is Conda due to its dependency on RDKit. Alongside RDKit, a number of different Python packages are required. An up-to-date list of package requirements are listed below as well as in the requirements.txt file:

numpy~=1.26.4

psycopg2~=2.8.6

pony~=0.7.14

lxml~=4.6.2

rdkit~=2024.09.3

wtforms~=3.2.1

flask~=2.3.3

werkzeug~=2.3.8

python-dateutil~=2.9.0.post0

requests~=2.32.3

scipy~=1.14.1

pillow~=11.0.0

pandas~=2.2.2

plotly~=5.24.1

chemcoord~=2.1.2

openbabel~=3.1.1

flask-wtf

flask-login



