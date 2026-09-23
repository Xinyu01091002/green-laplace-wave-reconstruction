"""Verify archived file integrity and runtime source boundaries, not physics."""
from pathlib import Path
import hashlib,json

root=Path(__file__).resolve().parent
manifest=json.loads((root/'package_files.json').read_text())
for item in manifest['files']:
    content=(root/item['path']).read_bytes()
    if item['normalization']=='LF':content=content.replace(b'\r\n',b'\n')
    assert hashlib.sha256(content).hexdigest()==item['sha256'],item['path']
gl=(root/'cpp/gl_eta44_fftw.cpp').read_text()
for forbidden in ['pair_stokes_corrections','third_order_stokes_correction',
                  'stokes_correction_diagonal','optimal_gain','shape_only']:
    assert forbidden not in gl,forbidden
assert 'GLRTI001' in gl
for file in root.glob('*.m'):
    body=file.read_text()
    assert 'C:/Users/' not in body and "'SWORD-VWA'" not in body,file.name
print(f"PASS: {len(manifest['files'])} package files and no-Stokes source boundary")
