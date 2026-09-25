"""Record an input-syntax-only correction while preserving prior generated fields."""
import hashlib
import json
import pathlib
import shutil
import sys

source, out = map(pathlib.Path, sys.argv[1:3])
assert source.is_dir() and not out.exists()
out.mkdir()
records = []
for phase in (0, 90, 180, 270):
    name = f'wavegroup_kpd1_akp012_phi{phase:03d}'
    folder = out / name
    folder.mkdir()
    shutil.copy2(source/name/'OceanWave3D.init', folder/'OceanWave3D.init')
    lines = (source/name/'OceanWave3D.inp').read_text().splitlines()
    assert lines[1] == '-1 0 1000' and len(lines[5].split()) == 1
    assert lines[12] == '0 0.0 0 X 0.0'
    lines[5] += ' 1000'
    selection = lines[9].split()
    selection[6] = '2'  # t=0 surface potential must be read from EP, not unsolved volume phi.
    lines[9] = ' '.join(selection)
    lines.insert(12, '0 0 0 0 0 0 0')
    (folder/'OceanWave3D.inp').write_text('\n'.join(lines)+'\n')
    hashes = {}
    for p in (source/name/'OceanWave3D.init', folder/'OceanWave3D.init', folder/'OceanWave3D.inp'):
        with p.open('rb') as data:
            hashes[str(p)] = hashlib.file_digest(data, 'sha256').hexdigest()
    assert hashes[str(source/name/'OceanWave3D.init')] == hashes[str(folder/'OceanWave3D.init')]
    records.append(dict(phase=phase, hashes=hashes))
manifest = dict(parent=str(source), initial_fields_mat=str(source/'initial_fields.mat'),
                changes=['Explicit rho=1000 on gravity line', 'Explicit disabled breaking-model line',
                         'Kinematics starts at time index 2; initial eta/psi comes from EP'],
                physical_fields_unchanged=True, cases=records)
(out/'repack.json').write_text(json.dumps(manifest, indent=2))
print(json.dumps(manifest, indent=2))
