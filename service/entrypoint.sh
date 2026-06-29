#!/bin/sh
# Container entrypoint: start the warm DG++Light inference server (loads the
# ~205 MB STRING Net-KNN index + integrators ONCE) in the background, then run
# the FastAPI app. Genome requests forward to the server via $DGPP_LIGHT_SERVER;
# while it is still warming up (socket not yet present) the sidecar falls back
# to in-process construction, so the API is never blocked on the model load and
# the server is a pure optimisation.
set -e

if [ -n "${DGPP_LIGHT_SERVER:-}" ]; then
  echo "[entrypoint] launching warm dgpp_light server -> ${DGPP_LIGHT_SERVER}"
  python3 /app/deepgo-plusplus/service/dgpp_server.py &
fi

exec uvicorn service.app:app --host 0.0.0.0 --port "${PORT:-8000}"
