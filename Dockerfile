FROM continuumio/miniconda3:4.9.2

RUN wget http://www.cibiv.at/software/pda/download/pda-1.0.3/pda-1.0.3-Linux.tar.gz \
    && tar -xf pda-1.0.3-Linux.tar.gz \
    && echo 'alias pda="/pda-1.0.3-Linux/bin/pda"' >> ~/.bashrc
RUN conda install -c bioconda cd-hit

COPY requirements.txt /temp/requirements.txt
RUN pip install -r /temp/requirements.txt

WORKDIR /down_sampling_analysis
COPY . .

RUN echo 'alias r4s="/down_sampling_analysis/docker_aux/Rate4Site/3.0/bin/rate4site_doublerep"' >> ~/.bashrc

ENV BASH_ENV=~/.bashrc
ENV PYTHONPATH=.

CMD ["python", "-m", "unittest", "discover", "./tests"]