# DeepGO-PlusPlus-Light webservice — strictly-CPU GO predictor.
# The image bundles the code, the frozen model JSONs and the DIAMOND binary.
# The large data assets (DIAMOND DB, net index, GO DAG/obo, train labels) are
# mounted at /assets at runtime (build them once with make_assets.sh).
FROM python:3.12-slim

# DIAMOND (CPU homology search) + curl to fetch the static binary
RUN apt-get update && apt-get install -y --no-install-recommends curl \
 && curl -L -o /tmp/diamond.tar.gz \
      https://github.com/bbuchfink/diamond/releases/download/v2.1.9/diamond-linux64.tar.gz \
 && tar -xzf /tmp/diamond.tar.gz -C /usr/local/bin diamond \
 && rm /tmp/diamond.tar.gz \
 && apt-get purge -y curl && apt-get autoremove -y && rm -rf /var/lib/apt/lists/*

WORKDIR /app
COPY service/requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

# code + frozen models (small, shipped in the repo)
COPY pipeline/apply_net_bridge.py /app/pipeline/apply_net_bridge.py
COPY service/predict.py service/app.py /app/service/
COPY models/deepgo_plusplus_light_fast.json \
     models/deepgo_plusplus_light_cpu.json /app/models/

ENV DGPP_ASSETS=/assets \
    DGPP_MODELS=/app/models \
    DGPP_THREADS=8
EXPOSE 8000
WORKDIR /app/service
# assets must be mounted at /assets (see make_assets.sh + README)
CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8000"]
