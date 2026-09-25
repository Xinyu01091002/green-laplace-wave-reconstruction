"""Freeze source, native binary and sixteen already-generated input pairs."""
import argparse
import datetime
import hashlib
import json
import pathlib
import re
import shutil
import subprocess
import tarfile

parser = argparse.ArgumentParser()
parser.add_argument('--project', required=True)
parser.add_argument('--inputs', required=True)
parser.add_argument('--run', required=True)
parser.add_argument('--source-archive', required=True)
parser.add_argument('--source-sha256', required=True)
parser.add_argument('--commit', required=True)
args = parser.parse_args()
project = pathlib.Path(args.project).resolve()
inputs = pathlib.Path(args.inputs).resolve()
root = pathlib.Path(args.run).resolve()
archive = pathlib.Path(args.source_archive).resolve()
assert root.parent == pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs')
assert not root.exists() and re.fullmatch(r'[a-zA-Z0-9_-]+', root.name)
def sha(path, payload=False):
    with pathlib.Path(path).open('rb') as data:
        if payload:
            data.readline(); data.readline()
        return hashlib.file_digest(data, 'sha256').hexdigest()
assert sha(archive) == args.source_sha256
with archive.open('rb') as data:
    archived_commit=subprocess.check_output(['git','get-tar-commit-id'],stdin=data,text=True).strip()
assert archived_commit == args.commit, 'Archive commit does not match the launch manifest'
matrix = json.loads((inputs/'case_matrix.json').read_text())
assert len(matrix['cases']) == 16
root.mkdir(parents=True)
source = root/'source'; source.mkdir()
with tarfile.open(archive) as bundle:
    bundle.extractall(source, filter='data')
shutil.copy2(archive, root/'source.tar')
source_hashes = {str(p.relative_to(source)):sha(p) for p in source.rglob('*') if p.is_file()}
(root/'bin').mkdir()
solver = root/'bin/ow3d'
shutil.copy2('/usr/local/bin/ow3d', solver)
assert sha(solver) == '36bd612d2655c26daa3a3146d7d3d291f0864d8cd3f23c3d2dd6992f7a26f614'
matlab = '/home/lxy/Desktop/matlabr2026a/bin/matlab'
assert pathlib.Path(matlab).is_file()
(root/'cases').mkdir();(root/'inputs').mkdir()
cases=[];payloads=set();family_hashes={}
for item in matrix['cases']:
    case = dict(item); identifier=case['case_id']
    assert re.fullmatch(r'[a-zA-Z0-9_]+', identifier)
    original=inputs/identifier; destination=root/'cases'/identifier; destination.mkdir()
    for name in ('OceanWave3D.inp','OceanWave3D.init'):
        shutil.copy2(original/name, destination/name)
    data_hash = sha(destination/'OceanWave3D.init', payload=True)
    assert data_hash not in payloads, 'Duplicate physical eta/psi payload; refusing redundant jobs'
    payloads.add(data_hash)
    case['input_hashes']={name:sha(destination/name) for name in ('OceanWave3D.inp','OceanWave3D.init')}
    case['initial_field_payload_sha256']=data_hash
    original_source=pathlib.Path(case['source_fields'])
    if not original_source.is_absolute():
        original_source=project/original_source
    family_dir=root/'inputs'/case['family']
    if not family_dir.exists():
        family_dir.mkdir()
        shutil.copy2(original_source, family_dir/'initial_fields.mat')
        family_hashes[case['family']]=sha(family_dir/'initial_fields.mat')
    case['source_fields']=str(family_dir/'initial_fields.mat')
    case['original_source_fields']=str(original_source)
    case['directory']=str(destination)
    (destination/'expected.json').write_text(json.dumps(case,indent=2))
    cases.append(case)
assert len(payloads)==16
with tarfile.open(root/'frozen-inputs.tar.gz','w:gz') as bundle:
    bundle.add(root/'inputs',arcname='inputs')
    for case in cases:
        for name in ('OceanWave3D.inp','OceanWave3D.init','expected.json'):
            bundle.add(root/'cases'/case['case_id']/name,arcname=f"cases/{case['case_id']}/{name}")
manifest=dict(run_id=root.name, created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
    source_commit=args.commit, source_archive_sha256=sha(root/'source.tar'), source_hashes=source_hashes,
    input_archive_sha256=sha(root/'frozen-inputs.tar.gz'), frozen_family_field_hashes=family_hashes,
    solver=str(solver),solver_sha256=sha(solver),matlab=matlab,cases=cases,concurrency=16,
    minimum_available_memory_bytes=560*1024**3,
    scientific_scope='kpd=1; same first-order modal amplitudes; one focused group and three random seeds, each four phases',
    boundary='Native OW3D walls; random results are not automatically equivalent to periodic HOS-Ocean',
    physical_accuracy_certified=False, automatic_fetch=False, automatic_cleanup=False,
    validation='Native completion, final time 220 s, full eta/phi strip integrity; no model-accuracy certification',
    runtime_timeout_seconds=None,
    mail='Existing remote run.sh recipient configuration; send once after all jobs and raw-output validations close')
(root/'manifest.json').write_text(json.dumps(manifest,indent=2))
(root/'source-commit.txt').write_text(args.commit+'\n')
(root/'status.txt').write_text('prepared\n')
(root/'resource-preflight.txt').write_text(pathlib.Path('/proc/meminfo').read_text())
environment=subprocess.run(['ldd',str(solver)],text=True,capture_output=True,check=True)
(root/'solver-libraries.txt').write_text(environment.stdout)
print(json.dumps(dict(run=str(root),case_count=len(cases),unique_payloads=len(payloads),
    source_archive_sha256=manifest['source_archive_sha256'],input_archive_sha256=manifest['input_archive_sha256']),indent=2))
