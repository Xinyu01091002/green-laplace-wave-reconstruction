"""Inventory a completed design attempt without promoting it to a validated run."""
import datetime
import hashlib
import json
import pathlib
import subprocess
import sys

base = pathlib.Path(sys.argv[1]).resolve()
repo = pathlib.Path(sys.argv[2]).resolve()
assert base.is_relative_to(repo/'results'/'ow3d_redesign')
destination = base/'MANIFEST.json'
assert not destination.exists(), 'Refusing to replace a frozen manifest.'
def read(relative):
    return json.loads((base/relative).read_text())
def digest(path):
    with path.open('rb') as data:
        return hashlib.file_digest(data, 'sha256').hexdigest()

paths = list((base/'source-v6').glob('*.m')) + list((base/'source-v6').glob('*.py'))
paths += list((base/'external_mf12').glob('*.m'))
paths += list((base/'wavegroup-v4').glob('*/OceanWave3D.*'))
paths += [base/'geometry-v1/first_order_design.mat', base/'wavegroup-v2/initial_fields.mat',
          base/'random-boundary-v1/random_boundary.mat', base/'wavegroup-v4/repack.json']
effective = read('wavegroup-v2/initialization.json')
effective['kinematics_physical_node_range'][6] = 2
effective['initial_output_source'] = 'EP_00000.bin at t=0; kinematics starts at dt=0.2 s'
effective['physical_fields_source'] = str(base/'wavegroup-v2/initial_fields.mat')
effective['effective_input_directory'] = str(base/'wavegroup-v4')
effective['propagation_validated'] = False
(base/'wavegroup-v4/effective_initialization.json').write_text(json.dumps(effective, indent=2))
manifest = dict(
    created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    repository_commit=subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip(),
    execution_base_commit=(base/'source-v6/BASE_COMMIT').read_text().strip(),
    source_snapshot='source-v6', solver='/usr/local/bin/ow3d',
    solver_sha256=digest(pathlib.Path('/usr/local/bin/ow3d')),
    hashes={str(p.relative_to(base)):digest(p) for p in paths},
    current_wavegroup_inputs='wavegroup-v4', physical_fields_source='wavegroup-v2/initial_fields.mat',
    geometry=read('geometry-v1/design.json'),
    native_output_audit=read('smoke-v4/initial-output-audit.json'),
    random_design=read('random-boundary-v1/random_boundary.json'),
    wavegroup_status='INPUTS_PREPARED; INITIAL_OUTPUT_AUDITED; TIME_ADVANCE_NOT_COMPLETED',
    random_status='UNIT_AMPLITUDE_FINITE_RANDOM_FIELD_PROTOTYPE; NORMALIZATION_AND_BOUNDARY_RUN_PENDING',
    production_started=False,
    observations=['First two-step OW3D attempt timed out after 300 s before first saved advance',
                  'Native t=0 kinematics phi differs from EP; production selection starts at index 2',
                  'No raw new fields fetched to the local workstation',
                  'Earlier parser failures and snapshots are preserved, not production inputs'])
destination.write_text(json.dumps(manifest, indent=2))
print(json.dumps({k:manifest[k] for k in ('repository_commit','wavegroup_status','random_status','production_started')},indent=2))
