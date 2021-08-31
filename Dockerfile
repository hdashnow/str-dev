FROM continuumio/miniconda3

#WORKDIR /app

# Create the environment
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda create -n strling-install strling

# Activate environment
RUN echo "conda activate strling-install" >> ~/.bashrc
SHELL ["/bin/bash", "--login", "-c"]

# Check install succeded and environment is active
RUN strling -h
RUN strling-outliers.py --help

