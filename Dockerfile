FROM julia:1.8.1

WORKDIR /app

ENV GKSwstype=100
ENV JULIA_NUM_THREADS=1
ENV MPLBACKEND=Agg
ENV OUTPUT_DIR=/app/output

RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt ./
RUN python3 -m pip install --no-cache-dir -r requirements.txt

COPY Project.toml Manifest.toml ./
RUN julia --project=/app -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'

COPY src ./src
COPY scripts ./scripts
COPY python ./python
COPY docker/entrypoint.sh /usr/local/bin/cml-model

RUN chmod +x /usr/local/bin/cml-model && mkdir -p /app/output

ENTRYPOINT ["cml-model"]
CMD ["figures"]
