"""Create a separately recorded two-step OW3D smoke case without changing its parents."""
import hashlib
import json
import pathlib
import shutil
import sys

source, out = map(pathlib.Path, sys.argv[1:3])
assert source.is_dir() and not out.exists()
out.mkdir()
shutil.copy2(source / 'OceanWave3D.init', out / 'OceanWave3D.init')
lines = (source / 'OceanWave3D.inp').read_text().splitlines()
assert lines[8].split()[1:] == ['20', '1', '1']
assert lines[9].split()[:6] == ['1', '1025', '1', '125', '133', '1']
original = {str(i): lines[i] for i in (4, 8, 9)}
lines[4] = '3 0.2 1 0.0 1 0.0'
lines[8] = '1 20 1 1'
lines[9] = '1 1025 1 125 133 1 1 3 1'
(out / 'OceanWave3D.inp').write_text('\n'.join(lines) + '\n')
files = [source/'OceanWave3D.init', source/'OceanWave3D.inp',
         out/'OceanWave3D.init', out/'OceanWave3D.inp']
hashes = {str(p): hashlib.file_digest(p.open('rb'), 'sha256').hexdigest() for p in files}
manifest = dict(purpose='Input/IO and resource smoke; not propagation validation',
                source=str(source), original_lines_zero_based=original,
                replaced_lines_zero_based={str(i): lines[i] for i in (4, 8, 9)},
                hashes=hashes, expected_final_time_s=.4)
(out/'preparation.json').write_text(json.dumps(manifest, indent=2))
print(json.dumps(manifest, indent=2))
