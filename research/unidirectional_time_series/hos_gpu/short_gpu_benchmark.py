"""Sequential GPU timing on the same two-second interval as the CPU benchmark."""
import argparse, ctypes, json, pathlib, subprocess, time
p=argparse.ArgumentParser();p.add_argument('--exe',required=True);p.add_argument('--out',required=True,type=pathlib.Path);a=p.parse_args()
a.out.mkdir(exist_ok=False)
nv=ctypes.CDLL('libnvidia-ml.so.1');assert nv.nvmlInit_v2()==0
h=ctypes.c_void_p();assert nv.nvmlDeviceGetHandleByIndex_v2(0,ctypes.byref(h))==0
class Memory(ctypes.Structure):
    _fields_=[('total',ctypes.c_ulonglong),('free',ctypes.c_ulonglong),('used',ctypes.c_ulonglong)]
def used():
    m=Memory();assert nv.nvmlDeviceGetMemoryInfo(h,ctypes.byref(m))==0;return m.used
results=[]
try:
    for family in ['low','high']:
        for precision,tolerance in [('double','1e-12'),('float','1e-9')]:
            name=family+'-'+precision;case=a.out/name;baseline=used();peak=baseline;start=time.monotonic()
            with (a.out/(name+'.log')).open('w') as log:
                proc=subprocess.Popen([a.exe,precision,'/root/hos-production-check-20260927T025257Z/'+family+'-adaptive.bin',str(case),'2',tolerance],stdout=log,stderr=subprocess.STDOUT)
                while proc.poll() is None:
                    peak=max(peak,used());time.sleep(.05)
            assert proc.returncode==0,(name,proc.returncode)
            row=json.loads((case/'summary.json').read_text());row.update(family=family,precision=precision,tolerance=float(tolerance),whole_process_seconds=time.monotonic()-start,baseline_vram_bytes=baseline,peak_vram_bytes=peak)
            results.append(row);(a.out/'summary.json').write_text(json.dumps(results,indent=2));print(json.dumps(row),flush=True)
finally:
    nv.nvmlShutdown()
