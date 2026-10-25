"""Public upload and results portal. Run one uvicorn process (one queue worker)."""
from contextlib import asynccontextmanager
import json
from pathlib import Path
import re
import shutil
import threading
import time
import uuid

from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.responses import FileResponse, JSONResponse
from service import jobs

STATIC = Path(__file__).parent / 'web'
PRESETS = {'bacteria', 'archaea', 'prokaryote', 'eukaryote', 'mag', 'metagenome', 'virus'}


@asynccontextmanager
async def lifespan(app):
    jobs.ROOT.mkdir(parents=True, exist_ok=True)
    jobs.STOP.clear()
    for path in jobs.ROOT.glob('*/status.json'):
        if jobs.state(path.parent)['status'] == 'running':
            jobs.state(path.parent, status='failed', stage='Interrupted by service restart',
                       error='The service restarted during this job. Please submit it again.')
    worker = threading.Thread(target=jobs.worker, daemon=True)
    worker.start()
    yield
    jobs.STOP.set()
    worker.join(timeout=5)


app = FastAPI(title='GSPA', lifespan=lifespan, docs_url=None, redoc_url=None)


@app.middleware('http')
async def headers(request, call_next):
    response = await call_next(request)
    response.headers['Referrer-Policy'] = 'no-referrer'
    response.headers['X-Content-Type-Options'] = 'nosniff'
    response.headers['Cache-Control'] = 'no-store'
    return response


@app.get('/')
def index():
    return FileResponse(STATIC / 'index.html')


@app.get('/health')
def health():
    ok = jobs.ready()
    return JSONResponse({'status': 'ready' if ok else 'unavailable', 'predictor': 'DeepGO++ Light',
                         'concurrency': 1, 'profile': 'full-v1', 'default_analyses': list(jobs.DEFAULT_ANALYSES)}, status_code=200 if ok else 503)


@app.post('/api/jobs', status_code=202)
def submit(file: UploadFile = File(...), kind: str = Form('gff3'), preset: str = Form('bacteria'),
           genome: UploadFile | None = File(None), disabled: str = Form(''), gram: str = Form('neg')):
    if kind not in ('gff3', 'genome', 'proteins') or preset not in PRESETS or gram not in ('neg', 'pos'):
        raise HTTPException(400, 'Select a valid input type and biological context.')
    if kind == 'genome' and preset == 'eukaryote':
        raise HTTPException(400, 'For eukaryotes, upload translated proteins; a splice-aware gene caller is required upstream.')
    disabled_set = set(filter(None, disabled.split(',')))
    if disabled_set - set(jobs.DEFAULT_ANALYSES):
        raise HTTPException(400, 'Unknown analysis option.')
    if not jobs.ready():
        raise HTTPException(503, 'Prediction assets are not ready. Please try again later.')
    with jobs.LOCK:
        active = sum(jobs.state(p.parent)['status'] in ('queued', 'running') for p in jobs.ROOT.glob('*/status.json'))
        if active >= 10:
            raise HTTPException(429, 'The queue is full. Please try again after a job finishes.')
        if shutil.disk_usage(jobs.ROOT).free < 8 * 1024**3:
            raise HTTPException(503, 'Storage is temporarily full. Please try again later.')
        job_id = uuid.uuid4().hex
        folder = jobs.ROOT / job_id
        folder.mkdir()
        try:
            def save(upload, name):
                size = 0
                with (folder / name).open('wb') as out:
                    while chunk := upload.file.read(1024 * 1024):
                        size += len(chunk)
                        if size > jobs.MAX_BYTES:
                            raise HTTPException(413, 'Each uncompressed input file must be at most 20 MB.')
                        out.write(chunk)
            name = {'genome': 'sample.fna', 'proteins': 'sample.faa', 'gff3': 'sample.gff3'}[kind]
            save(file, name)
            if genome is not None and genome.filename:
                if kind != 'gff3':
                    raise ValueError('A companion genome FASTA is only used with GFF3 input.')
                save(genome, 'sample.fna')
            if kind == 'gff3':
                from service.inputs import gff3
                count, length = gff3(folder)
            else:
                count, length = jobs.validate(folder / name, kind)
            data = jobs.state(folder, id=job_id, status='queued', stage='Waiting for a worker',
                created=time.time(), kind=kind, preset=preset, sequences=count, residues=length,
                analyses={key: key not in disabled_set for key in jobs.DEFAULT_ANALYSES}, gram=gram,
                filename=Path(file.filename or 'FASTA').name[:180])
        except Exception as exc:
            shutil.rmtree(folder)
            if isinstance(exc, ValueError):
                raise HTTPException(400, str(exc))
            raise
    return data


def folder_for(job_id):
    if not re.fullmatch('[0-9a-f]{32}', job_id):
        raise HTTPException(404, 'Job not found.')
    folder = jobs.ROOT / job_id
    if not (folder / 'status.json').is_file():
        raise HTTPException(404, 'Job not found or its seven-day retention period has ended.')
    return folder


@app.get('/api/jobs/{job_id}')
def status(job_id: str):
    return jobs.state(folder_for(job_id))


@app.get('/api/jobs/{job_id}/{artifact}')
def artifact(job_id: str, artifact: str):
    folder = folder_for(job_id)
    paths = {'browser': folder / 'gspa_out/sample_browser.html',
             'download': folder / 'results.zip', 'log': folder / 'run.log',
             'model': folder / 'gspa_out/metabolism/metabolic_model.xml',
             'fluxes': folder / 'gspa_out/metabolism/reactions.tsv',
             'analyses': folder / 'gspa_out/analyses.json',
             'annotation': folder / 'sample.gff3', 'sequence': folder / 'sample.fna'}
    if artifact not in paths or not paths[artifact].is_file():
        raise HTTPException(404, 'This result is not available yet.')
    if artifact != 'log' and jobs.state(folder)['status'] != 'complete':
        raise HTTPException(409, 'The job has not completed.')
    return FileResponse(paths[artifact], filename=None if artifact == 'browser' else paths[artifact].name,
                        media_type='text/html' if artifact == 'browser' else None)


@app.get('/examples/{name}')
def example(name: str):
    return example_artifact(name, 'browser')


@app.get('/examples/{name}/{kind}')
def example_artifact(name: str, kind: str):
    path = jobs.ROOT / 'examples.json'
    examples = json.loads(path.read_text()) if path.exists() else {}
    if name not in examples:
        raise HTTPException(404, 'Example not available.')
    return artifact(examples[name], kind)
