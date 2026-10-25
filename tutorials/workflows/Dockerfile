# A small, CPU-capable teaching environment built from an unconfigured OS image.
FROM ubuntu:24.04 AS build
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates openjdk-21-jdk-headless python3 python3-venv \
    diamond-aligner prodigal curl unzip \
    && rm -rf /var/lib/apt/lists/*
WORKDIR /src
COPY . .
RUN --mount=type=cache,target=/root/.gradle,sharing=locked ./gradlew :gspa-cli:shadowJar --no-daemon --max-workers=2 \
    && mkdir /opt/gspa && cp gspa-cli/build/libs/gspa-2.1.0.jar /opt/gspa/gspa-2.1.0.jar \
    && cp gspa /opt/gspa/gspa
ENV PATH="/opt/gspa:${PATH}"
ENV GSPA_JVM_OPTS="-Xmx4g -XX:ActiveProcessorCount=2"
WORKDIR /work
CMD ["gspa", "--help"]
