# extraCellularRNA

find DESeq data sets

anything with de.seq or de-seq is going to be de seq output of either the normalized counts or differential expression variety
 ```
 find /public/groups/kimlab  -name "*de.seq*" -print |grep -v "permission denied"
 ```
 
conda references
- [https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html)
- [cheat sheet](https://docs.conda.io/projects/conda/en/4.6.0/_downloads/52a95608c49671267e40c689e0bc00ca/conda-cheatsheet.pdf)


start the env

```
cd extraCellularRNA
conda activate extraCellularRNA
export PYTHONPATH="${PYTHONPATH}:`pwd`/src"
```

create the extraCellularRNA environment from yaml file

```
conda env create -f environment.yml
pip install tensorflow
```

updating dependencies
see [exporting-an-environment-file-across-platforms](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#exporting-an-environment-file-across-platforms)

```
cd extraCellularRNA
conda activate extraCellularRNA
conda env export --from-history > environment.yml
```

Running Unit test

```
cd extraCellularRNA
conda activate extraCellularRNA
export PYTHONPATH="${PYTHONPATH}:`pwd`/src"
cd src/test
python -m unittest discover .
```
