# init a base image
FROM continuumio/miniconda3
SHELL ["/bin/bash", "-c"]
RUN apt-get update
RUN apt-get install -y openbabel \
    && apt-get -y install libpq-dev gcc \
    && pip install psycopg2
# define the present working directory
WORKDIR /app
# copy the contents into the working dir
ADD . /app
COPY environment.yml .
RUN conda env create -f environment.yml
# install MOPAC
RUN conda install conda-forge::mopac
RUN chmod 777 /app/mopac/
# RUN chmod +x /app/mopac/MOPAC2016.exe
ENV LD_LIBRARY_PATH=/app/mopac
# Create the environment:

# Make RUN commands use the new environment:
SHELL ["conda", "run", "-n", "CH_deploy", "/bin/bash", "-c"]
# install RDKit
RUN conda install -c conda-forge rdkit
# install chemcoord
RUN conda install -c conda-forge pymatgen
# run pip to install the dependencies of the flask app
RUN pip install -r requirements.txt
COPY requirements.txt requirements.txt

# The code to run when container is started:
COPY . .
EXPOSE 5000
# ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "CH_deploy", "python", "app.py"]
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "CH_deploy", "python", "-m", "flask", "run"]