# GSPA - Genome-Scale Protein Annotation
# Multi-stage build: build JAR then create runtime image with bioinformatics tools

# --- Stage 1: Build GSPA JAR ---
FROM eclipse-temurin:21-jdk AS builder

WORKDIR /build
COPY gradlew settings.gradle.kts build.gradle.kts gradle.properties ./
COPY gradle/ gradle/
COPY gspa-core/ gspa-core/
COPY gspa-predictors/ gspa-predictors/
COPY gspa-cli/ gspa-cli/

RUN ./gradlew :gspa-cli:installDist --no-daemon

# --- Stage 2: Build gapsmith (Rust) ---
# Pinned version; bump when a new gapsmith release with compatible
# *-Reactions.tbl / *-Pathways.tbl output is needed.
FROM rust:1.82-bookworm AS gapsmith-builder
ARG GAPSMITH_REF=v0.1.0
RUN apt-get update && apt-get install -y --no-install-recommends \
        cmake clang libclang-dev pkg-config git \
    && rm -rf /var/lib/apt/lists/*
RUN git clone --depth 1 --branch ${GAPSMITH_REF} \
        https://github.com/bio-ontology-research-group/gapsmith.git /src/gapsmith
WORKDIR /src/gapsmith
RUN cargo build --release --bin gapsmith

# --- Stage 3: Runtime ---
FROM eclipse-temurin:21-jre

LABEL maintainer="Robert Hoehndorf <robert.hoehndorf@kaust.edu.sa>"
LABEL description="GSPA - Genome-Scale Protein Annotation pipeline"

# Install bioinformatics tools.
# Note: minced is intentionally not bundled — it was removed from
# Ubuntu 24.04 (Noble) repos. CRISPR detection is wired through the
# dedicated biocontainer in gspa-nf (`quay.io/biocontainers/minced`),
# so this orchestration image doesn't need it baked in.
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git \
    python3 python3-pip \
    prodigal \
    hmmer \
    ncbi-blast+ \
    barrnap \
    && rm -rf /var/lib/apt/lists/*

# Install DIAMOND
RUN wget -q https://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz \
    && tar xzf diamond-linux64.tar.gz && mv diamond /usr/local/bin/ && rm diamond-linux64.tar.gz

# Install MMseqs2
RUN wget -q https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz \
    && tar xzf mmseqs-linux-avx2.tar.gz && mv mmseqs/bin/mmseqs /usr/local/bin/ \
    && rm -rf mmseqs mmseqs-linux-avx2.tar.gz

# Install pyrodigal (faster Prodigal)
RUN pip3 install --break-system-packages pyrodigal

# Install Python tools for optional FBA
RUN pip3 install --break-system-packages cobra || true

# Copy GSPA distribution
COPY --from=builder /build/gspa-cli/build/install/gspa-cli /opt/gspa
ENV PATH="/opt/gspa/bin:${PATH}"

# Install gapsmith (Phase 11 default pathway tool) + vendored curation data.
# Large public tables (SEED, MNXref) and sequence DB (~2 GB) are NOT baked
# in — fetch them at run time with `gapsmith update-data` / `update-sequences`,
# then point `--data-dir` at the result.
COPY --from=gapsmith-builder /src/gapsmith/target/release/gapsmith /usr/local/bin/gapsmith
COPY --from=gapsmith-builder /src/gapsmith/data /opt/gapsmith/data

# Default working directory
WORKDIR /data

ENTRYPOINT ["gspa-cli"]
CMD ["--help"]
