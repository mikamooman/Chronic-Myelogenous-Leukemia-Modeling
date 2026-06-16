# CML Model Reproduction Code

This repository contains Julia code for reproducing stochastic simulations, model fitting, and figures for a model of chronic myeloid leukemia (CML) initiation.

The recommended way to run the code is with Docker. Docker installs the required Julia environment and package versions inside a container, so users do not need to manually configure Julia packages.

## Requirements

Install Docker Desktop:

* macOS: Docker Desktop for Mac
* Windows: Docker Desktop with WSL 2 backend
* Linux: Docker Engine or Docker Desktop

Check that Docker is working:

```bash
docker --version
docker info
```

## Build the Docker image

From the top-level repository folder:

```bash
docker build -t cml-model:paper .
```

## Reproduce all figures

```bash
docker run --rm \
  --env JULIA_NUM_THREADS=1 \
  -v "$(pwd)/output:/app/output" \
  cml-model:paper figures
```

Generated files will be written to the local `output/` folder.

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

The full fit may take substantially longer than the demonstration fit.

## Repository structure

```text
.
├── Dockerfile
├── Project.toml
├── Manifest.toml
├── README.md
├── src/
├── scripts/
├── docker/
└── output/
```

## Notes on threading

The Docker commands above use:

```bash
--env JULIA_NUM_THREADS=1
```

This is recommended for reproducibility because the stochastic simulations use random sampling.

## Running without Docker

Users with Julia installed can run:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
julia --project=. scripts/reproduce_all_figures.jl
```

Individual scripts can also be run directly:

```bash
julia --project=. scripts/reproduce_survival_figure.jl
julia --project=. scripts/reproduce_chimerism_figure.jl
julia --project=. scripts/demo_fit.jl
julia --project=. scripts/fit_model.jl
```

## Citation

If this repository is used, please cite the associated manuscript.

```text
[Add citation here]
```

## License

```text
[Add license here]
```
