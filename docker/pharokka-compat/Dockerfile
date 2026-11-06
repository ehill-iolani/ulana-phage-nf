# Local-only Pharokka image for Apple-silicon Macs. The stock biocontainer
# segfaults under Docker's amd64 emulation: polars 1.42.1's default runtime
# (polars-runtime-32) uses AVX2, which the emulator doesn't provide. The
# "compat" runtime is the same polars without it. Real x86 machines (Google
# Batch, CI) use the stock image; this one is never mirrored or pushed.
#
#   docker build --platform linux/amd64 -t pharokka-local:1.10.1 docker/pharokka-compat
#   nextflow run main.nf -profile test,docker --pharokka_container pharokka-local:1.10.1
FROM --platform=linux/amd64 quay.io/biocontainers/pharokka:1.10.1--pyhdfd78af_0
RUN pip install --no-cache-dir polars-runtime-compat==1.42.1 && pip uninstall -y polars-runtime-32
