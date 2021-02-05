FROM continuumio/miniconda3:4.9.2 AS base

# install depndency programs
RUN wget http://www.cibiv.at/software/pda/download/pda-1.0.3/pda-1.0.3-Linux.tar.gz \
    && tar -xf pda-1.0.3-Linux.tar.gz
RUN conda install -c bioconda cd-hit
RUN conda install -c bioconda mafft
RUN conda install -c bioconda prank
RUN conda install -c bioconda raxml

# copy code and toy data
COPY /docker_aux/rate4site_doublerep rate4site_doublerep
RUN ["chmod", "777", "rate4site_doublerep"]

# set aliases
RUN echo 'alias pda="~/pda-1.0.3-Linux/bin/pda"' >> ~/.bashrc
RUN echo 'alias rate4site="~/rate4site_doublerep"' >> ~/.bashrc

# install python packages
COPY requirements.txt /temp/requirements.txt
RUN pip install -r /temp/requirements.txt

ENV BASH_ENV=~/.bashrc

WORKDIR /src
COPY /src .
###########START NEW IMAGE : DEBUGGER ###################
FROM base AS debug

COPY /data/test /data/test
CMD ["python main.py"]

###########START NEW IMAGE: PRODUCTION ###################

FROM base AS prod

ENTRYPOINT ["python", "main.py"]
CMD ["--input_path=input.json"]