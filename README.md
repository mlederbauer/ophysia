<p align="center">
  <img src="assets/banner.png" width="200">
</p>

[![Code style: ruff-format](https://img.shields.io/badge/code%20style-ruff_format-6340ac.svg)](https://github.com/astral-sh/ruff)
![Coverage Status](https://raw.githubusercontent.com/mlederbauer/ophysia/main/coverage-badge.svg)

<h1 align="center">
ophysia
</h1>

<br>


ORCA 6.0 and related scripts

## 🔥 Usage

The goal of this project is to provide simple and easy-to-use scripts for generating ORCA input and scraping ORCA output files. Made to work with the ETH Euler Supercomputing cluster, can be extended to any other cluster (and local).

```
(ophysia) $ ophysia create \
    --smiles "CCO" \
    --charge 0 \
    --multiplicity 1 \
    --basis "def2-SVP" \
    --functional "PBE0" \
    --jobtype "opt" \
    --memory 4 \
    --cores 4 \
    --time 24
```

## 👩‍💻 Installation

Create a new environment, you may also give the environment a different name. 

```
conda create -n ophysia python=3.10 
```

```
conda activate ophysia
```

If you need jupyter lab, install it 

```
(ophysia) $ pip install jupyterlab
```


## 🛠️ Development installation

To install, run

```
(ophysia) $ pip install -e ".[test,doc]"
```

When adding new packages to the file, add them to the `pyproject.toml` file:

```
# pyproject.toml
...
[project]
name = "ophysia"
...
dependencies = [
    "foo-package==1.0.1",
    "foo-package-from-repo @ git+https://github.com/foo-username/foo-package-from-repo.git@main"
]
```

To run style checks:

```
(ophysia) $ pip install pre-commit
(ophysia) $ pre-commit run -a
```

### Run style checks, coverage, and tests

```
(ophysia) $ pip install tox
(ophysia) $ tox
```

### Generate coverage badge

Works after running `tox`

```
(ophysia) $ pip install "genbadge[coverage]"
(ophysia) $ genbadge coverage -i coverage.xml
```


