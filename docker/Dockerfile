# GSPA - Genome-Scale Protein Annotation
# Multi-stage build: build JAR then create runtime image with bioinformatics tools

# --- Stage 1: Build ---
FROM eclipse-temurin:21-jdk AS builder

WORKDIR /build
COPY gradlew settings.gradle.kts build.gradle.kts gradle.properties ./
COPY gradle/ gradle/
COPY gspa-core/ gspa-core/
COPY gspa-predictors/ gspa-predictors/
COPY gspa-cli/ gspa-cli/

RUN ./gradlew :gspa-cli:installDist --no-daemon

# --- Stage 2: Runtime ---
FROM eclipse-temurin:21-jre

LABEL maintainer="Robert Hoehndorf <robert.hoehndorf@kaust.edu.sa>"
LABEL description="GSPA - Genome-Scale Protein Annotation pipeline"

# Install bioinformatics tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    wget curl git \
    python3 python3-pip \
    prodigal \
    hmmer \
    ncbi-blast+ \
    minced \
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

# Default working directory
WORKDIR /data

ENTRYPOINT ["gspa-cli"]
CMD ["--help"]
