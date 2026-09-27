"""Report an observed throughput window, not a matched-work benchmark."""
import json,pathlib
r=pathlib.Path(__file__).resolve().parent
rows=[]
for line in (r/'parallel-samples.jsonl').read_text().splitlines():
    try:a=json.loads(line)
    except json.JSONDecodeError:continue
    if len(a['phases'])==4 and all('time' in p for p in a['phases'].values()):rows.append(a)
a,b=rows[0],rows[-1];elapsed=b['monotonic_seconds']-a['monotonic_seconds'];assert elapsed>=60
rates={n:(b['phases'][n]['time']-a['phases'][n]['time'])/elapsed for n in a['phases']}
old=json.loads((r/'parallel-adoption.json').read_text())['solver_progress'];baseline=old['time']/old['wall_seconds']
result={'window_wall_seconds':elapsed,'per_phase_simulated_seconds_per_wall_second':rates,'aggregate_simulated_seconds_per_wall_second':sum(rates.values()),'prior_single_phase_average':baseline,'observed_aggregate_ratio':sum(rates.values())/baseline,'device_used_bytes':b['device_used_bytes'],'last_progress':b['phases'],'limitation':'Before/after observation at different phase/time positions; not a matched-work serial-versus-parallel benchmark.'}
(r/'parallel-observation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
