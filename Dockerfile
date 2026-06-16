FROM julia:1.8.1

WORKDIR /app

ENV GKSwstype=100
ENV JULIA_NUM_THREADS=1

COPY Project.toml Manifest.toml ./

RUN julia --project=/app -e \
    'using Pkg; Pkg.instantiate(); Pkg.precompile()'

COPY src ./src
COPY scripts ./scripts
COPY docker/entrypoint.sh /usr/local/bin/cml-model

RUN chmod +x /usr/local/bin/cml-model \
    && mkdir -p /app/output

ENTRYPOINT ["cml-model"]
CMD ["figures"]
