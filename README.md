# CML Model Reproduction Code

This repository contains the Julia and Python code used to reproduce simulations, model fitting, and figures for a mathematical model of chronic myeloid leukemia (CML) initiation. The model describes the stochastic dynamics of normal and leukemic hematopoietic cell populations and examines how feedback regulation and initial leukemic chimerism influence disease progression.

Docker is the recommended way to run the repository. It provides the required Julia and Python environments, package versions, and system dependencies.


## Requirements

Install Docker Desktop: https://docs.docker.com/get-started/get-docker/

- windows will also need to install WSL: https://learn.microsoft.com/en-us/windows/wsl/install

Check that Docker is working:

```bash
docker --version
docker info
```

## Tested

so far tested on macOS Tahoe 26.5.1


##Download from git

```bash
git clone https://github.com/mikamooman/Chronic-Myelogenous-Leukemia-Modeling.git
```

```bash
cd Chronic-Myelogenous-Leukemia-Modeling
```

## Build the Docker image

From the top-level repository folder:

```bash
docker build -t cml-model:paper .
```
this should take ~ 5-10 mins


## Reproduce CML free survival plots and Chimerism vs time plots

if using bash try:

```bash
docker run --rm \
  --env JULIA_NUM_THREADS=1 \
  -v "$(pwd)/output:/app/output" \
  cml-model:paper figures
```

if using powershell try:
```bash
docker run --rm --env JULIA_NUM_THREADS=1 -v "${pwd}/output:/app/output" cml-model:paper figures
```


Generated files will be written to the local `output/` folder.

expect this to take ~10-30 mins

## Run a small demo fit and plot the results

```bash
docker run --rm \
  --env JULIA_NUM_THREADS=1 \
  -v "$(pwd)/output:/app/output" \
  cml-model:paper demo-fit
```

if using powershell try:
```bash
docker run --rm --env JULIA_NUM_THREADS=1 -v "${pwd}/output:/app/output" cml-model:paper demo-fit
```

Generated files will be written to the local `output/` folder.

expect this to take ~ 30-45 mins

## Recreate supplementary CHIP figures

```bash
docker run --rm \
  --env JULIA_NUM_THREADS=1 \
  -v "$(pwd)/output:/app/output" \
  cml-model:paper supplemental-python

```

if using powershell try:
```bash
docker run --rm --env JULIA_NUM_THREADS=1 -v "${pwd}/output:/app/output" cml-model:paper supplemental-python
```


this shouldnt take long < 5 minutes

## Available workflows

Reproduce all figures:

```bash
docker run --rm \
  --env JULIA_NUM_THREADS=1 \
  -v "$(pwd)/output:/app/output" \
  cml-model:paper figures
```

Reproduce only the CML-free survival figure:

```bash
docker run --rm \
  --env JULIA_NUM_THREADS=1 \
  -v "$(pwd)/output:/app/output" \
  cml-model:paper survival
```

Reproduce only the chimerism figure:

```bash
docker run --rm \
  --env JULIA_NUM_THREADS=1 \
  -v "$(pwd)/output:/app/output" \
  cml-model:paper chimerism
```

Run a short demonstration fit:

```bash
docker run --rm \
  --env JULIA_NUM_THREADS=1 \
  -v "$(pwd)/output:/app/output" \
  cml-model:paper demo-fit
```

Run the full model fit:

```bash
docker run --rm \
  --env JULIA_NUM_THREADS=1 \
  -v "$(pwd)/output:/app/output" \
  cml-model:paper full-fit
```

Note: full model fit will take hours possibly a day or two.

Reproduce the supplementary CHIP figures

```bash
docker run --rm \
  --env JULIA_NUM_THREADS=1 \
  -v "$(pwd)/output:/app/output" \
  cml-model:paper supplemental-python
```

The full fit may take substantially longer than the demonstration fit.

## Notes on threading

The Docker commands above use:

```bash
--env JULIA_NUM_THREADS=1
```

Multi - threading is not currently supported

## Language and package versions

julia: 1.8.1
Optim: 1.11.0
Distributions: 0.25.118
Plots: 1.40.20
StatsBase: 0.34.4

Python: 3.12.3
Numpy: 2.3.4
Scipy: 1.16.3
Matplotlib: 3.9.2.



```

## Citation

If this repository is used, please cite the associated manuscript.

```text
[Add citation here]
```

## License

```text

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

```
