"""Summarize allocation *samples*, not retained-object sizes, from JFR JSON."""
import collections, hashlib, json, sys
from pathlib import Path
import ijson
source=Path(sys.argv[1])
counts=collections.Counter(); classes=collections.Counter(); sites=collections.Counter(); heaps=[]; writes=[]
with source.open('rb') as stream:
    for event in ijson.items(stream,'recording.events.item',use_float=True):
        kind=event['type']; v=event['values']; counts[kind]+=1
        if kind=='jdk.ObjectAllocationSample':
            weight=v['weight']; classes[v['objectClass']['name']]+=weight
            frames=(v.get('stackTrace') or {}).get('frames',[])
            frame=next((f for f in frames if f['method']['type']['name'].startswith('gspa/')),None)
            key='no GSPA frame in sampled stack' if frame is None else f"{frame['method']['type']['name']}.{frame['method']['name']}:{frame['lineNumber']}"
            sites[key]+=weight
        elif kind=='jdk.GCHeapSummary':
            heaps.append({k:v[k] for k in ('startTime','gcId','when','heapUsed')})
        elif kind=='jdk.FileWrite':
            writes.append({k:v.get(k) for k in ('startTime','duration','path','bytesWritten')})
with source.open('rb') as stream:
    digest=hashlib.file_digest(stream,'sha256').hexdigest()
print(json.dumps(dict(source_sha256=digest, event_counts=counts,
    interpretation='Allocation weights estimate allocated bytes over time, including dead temporary objects; they do not measure retained heap. Before/after GC heap values are observations, not a live-object census.',
    sampled_allocations_by_class=classes.most_common(20),sampled_allocations_by_first_gspa_frame=sites.most_common(25),
    max_observed_heap_bytes=max((h['heapUsed'] for h in heaps),default=None),
    gc_heap_summaries=heaps,file_writes=writes),indent=2))
