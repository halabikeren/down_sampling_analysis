FROM continuumio/miniconda3:4.9.2 AS base

# install depndency programs
RUN wget http://www.cibiv.at/software/pda/download/pda-1.0.3/pda-1.0.3-Linux.tar.gz \
    && tar -xf pda-1.0.3-Linux.tar.gz
RUN conda install -c bioconda cd-hit
RUN conda install -c bioconda mafft
RUN conda install -c bioconda prank
RUN conda install -c bioconda raxml

# set aliases
RUN echo 'alias pda="/pda-1.0.3-Linux/bin/pda"' >> ~/.bashrc
RUN echo 'alias rate4site="/down_sampling_analysis/docker_aux/Rate4Site/3.0/bin/rate4site_doublerep"' >> ~/.bashrc

# install python packages
COPY requirements.txt /temp/requirements.txt
RUN pip install -r /temp/requirements.txt

# copy code and toy data
WORKDIR /app/down_sampling_analysis
COPY . .
RUN chmod 777 /down_sampling_analysis/docker_aux/rate4site_doublerep

ENV BASH_ENV=~/.bashrc
ENV PYTHONPATH=.

###########START NEW IMAGE : DEBUGGER ###################
FROM base AS debug

CMD ["python", "-m", "unittest", "discover", "./tests"]

###########START NEW IMAGE: PRODUCTION ###################
FROM base AS prod

RUN mkdir /down_sampling_analysis/workdir

