![unit tests](https://github.com/halabikeren/down_sampling_analysis/workflows/unit%20tests/badge.svg)

[![https://img.shields.io/docker/pulls/halabikeren/sap.svg](https://img.shields.io/docker/pulls/halabikeren/dsap.svg)](https://hub.docker.com/repository/docker/halabikeren/dsap)
 
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/5160)

This project compares the performance of several down-samplnig methods in aim of implementing a single optimal one

# Input

A folder with a json file and a sequence data file. See example [here](https://github.com/halabikeren/down_sampling_analysis/tree/master/data/test). For additional documentation on the available input parameters, go [here](https://github.com/halabikeren/down_sampling_analysis/blob/master/data/input_json_parameters.txt).

# Docker usage
```
docker pull halabikeren/down_sampling_analysis_prod:latest
docker run -v <path_to_input_filder>:/data/ dsap --input_path=/data/<json_filename>
```


# Singularity usage
```
singularity pull --name dsap.sif shub://halabikeren/down_sampling_analysis
singularity run --bind <path_to_input_filder>:/data/ --writable-tmpfs dsap.sif --input_path=/data/<json_filename>
```

# Direct usage
```
python /src/main.py --input_path=<path_to_input_filder><json_filename>
```


