#!/usr/bin/env bash

set -euo pipefail

COMMAND="${1:-figures}"

case "$COMMAND" in
    figures)
        exec julia --project=/app /app/scripts/reproduce_all_figures.jl
        ;;

    survival)
        exec julia --project=/app /app/scripts/reproduce_survival_figure.jl
        ;;

    chimerism)
        exec julia --project=/app /app/scripts/reproduce_chimerism_figure.jl
        ;;

    demo-fit)
        exec julia --project=/app /app/scripts/demo_fit.jl
        ;;

    full-fit)
        exec julia --project=/app /app/scripts/fit_model.jl
        ;;

    bash)
        exec /bin/bash
        ;;


  supplemental-python)
    exec python3 /app/python/supplemental_chip_figures.py
    ;;

    *)
        echo "Unknown command: $COMMAND"
        echo "Available commands:"
        echo "  figures"
        echo "  survival"
        echo "  chimerism"
        echo "  demo-fit"
        echo "  full-fit"
        echo "  bash"
        exit 1
        ;;
esac
