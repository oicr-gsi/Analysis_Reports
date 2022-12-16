# Analysis_Reports

## Introduction ##
This repository provides a script `ar.py` to generate analysis reports to accompany released data.


## Getting Started ##

### Installation ###

Clone this repository into a folder on the cluster.
```
git clone git@github.com:oicr-gsi/Analysis_Reports.git
```

### Dependenices ###

The file `requirements.txt` has the dependenices required for the script to work. It is recomended that you install these in a new python virtual environment. One option is to use [Anaconda](https://www.anaconda.com/), a popular package management software that was used in the development of this repository. 

After installing Anaconda and setting up a new virtual environment, run the following command in the root directory of this repo
```
conda install --file requirements.txt
```

If you opt to use a different package manager like `Pyenv`, you will most likely run 
```
pip install -r requirements.txt
```

## Usage ##

Two examples of using the script are as follows

1. ```
python3 generate.py -i infile.json -o outfile.pdf
```
This will pull data from `infile.json` and, using data from production, will create an Analysis Report in the current directory named `outfile.pdf`

2. ```
python3 generate.py -i infile.json -o /reports/outfile.pdf --stage
```
This will pull data from `infile.json` and, using data from stage, will create an Analysis Report in the folder `reports` directory named `outfile.pdf`.


###Parameters###

| argument | abbreviation| purpose | required/optional | default                                   |
| ------- | ------- | ------------------------------------------ |
| --infile | -i| The input file to be read| required | ar_input.json |
| --outfile | -o | Name of the report generated | required | "ar_report.pdf"  |
| --staging, --stage | |If used, data will be pulled from stage| optional | Leaving out the flag will pull data from production |



#### WIP Basic json structure ####

The basic structure is organized with sample ids, and workflow names and ids.

```
{"cases":
   {workflow_1:
      {"workflow_id": "19168526", "workflow_version":"2.0.2"},
    workflow_2:
      {"workflow_id":"16962244", "workflow_version":"2.0.2"}
   }
}
```

This will create the following directory structure under the directory specified by the `-p` parameter:

```
sample_id
 |- workflow_1
    |- file 1
    |- file 2
    |- ...
 |- workflow_2
    |- file 1
    |- file 2
    |- ...
```
