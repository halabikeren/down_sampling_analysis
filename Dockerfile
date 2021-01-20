FROM continuumio/miniconda3:4.9.2

RUN wget http://www.cibiv.at/software/pda/download/pda-1.0.3/pda-1.0.3-Linux.tar.gz \
    && tar -xf pda-1.0.3-Linux.tar.gz \
    && alias pda="/pda-1.0.3-Linux/bin/pda"
RUN conda install -c bioconda cd-hit

COPY requirements.txt /temp/requirements.txt
RUN pip install -r /temp/requirements.txt

WORKDIR /app
COPY . .
CMD ["python"]