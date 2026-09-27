"""Remote sequential measurement orchestration; never changes GPU settings."""
import ctypes
import json
import os
import pathlib
import subprocess
import time

root=pathlib.Path(__file__).resolve().parent
out=root/'measurements';out.mkdir(exist_ok=False)
nvml=ctypes.CDLL('libnvidia-ml.so.1')
assert nvml.nvmlInit_v2()==0
handle=ctypes.c_void_p();assert nvml.nvmlDeviceGetHandleByIndex_v2(0,ctypes.byref(handle))==0
class Memory(ctypes.Structure):
    _fields_=[('total',ctypes.c_ulonglong),('free',ctypes.c_ulonglong),('used',ctypes.c_ulonglong)]
def memory():
    m=Memory();assert nvml.nvmlDeviceGetMemoryInfo(handle,ctypes.byref(m))==0;return m.used
def temperature():
    t=ctypes.c_uint();assert nvml.nvmlDeviceGetTemperature(handle,0,ctypes.byref(t))==0;return t.value
metadata={}
for name,cmd in [('cpu',['lscpu']),('gpu',['nvidia-smi','-q']),('compiler',['g++','--version'])]:
    r=subprocess.run(cmd,text=True,capture_output=True);metadata[name]=r.stdout+r.stderr
(out/'hardware.json').write_text(json.dumps(metadata,indent=2))
env=dict(os.environ,OMP_NUM_THREADS='8',OMP_DYNAMIC='FALSE',OMP_PROC_BIND='close',OMP_PLACES='cores')
rows=[]
for nx,ny in [(512,256),(2048,1024)]:
    for precision in ['float','double']:
        for backend in ['cpu','gpu']:
            name=f'{backend}-{precision}-{nx}x{ny}'
            command=['taskset','-c','0-7',str(root/'build'/f'hos_benchmark_{backend}'),precision,str(nx),str(ny),'8','3',str(out/(name+'.bin'))]
            baseline=memory();peak=baseline;temp0=temperature();temp_peak=temp0
            print('Starting '+name,flush=True);start=time.monotonic();last=start
            with (out/(name+'.log')).open('w') as log:
                p=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env)
                while p.poll() is None:
                    peak=max(peak,memory());temp_peak=max(temp_peak,temperature())
                    if time.monotonic()-last>20:
                        print(name+' running; elapsed %.1f s'%(time.monotonic()-start),flush=True);last=time.monotonic()
                    if time.monotonic()-start>1200:
                        p.terminate();p.wait();raise TimeoutError(name)
                    time.sleep(.05)
            if p.returncode:
                print((out/(name+'.log')).read_text(),flush=True);raise RuntimeError(name+' failed')
            lines=(out/(name+'.log')).read_text().splitlines()
            row=json.loads(next(x for x in reversed(lines) if x.startswith('{')))
            row.update(name=name,nvml_baseline_bytes=baseline,nvml_peak_bytes=peak,
                       nvml_increment_peak_bytes=peak-baseline,temp_before=temp0,temp_peak_sampled=temp_peak,
                       temp_after=temperature(),nvml_sample_interval_seconds=.05)
            rows.append(row);(out/'timings.json').write_text(json.dumps(rows,indent=2))
            print(json.dumps(row),flush=True)
checks=[]
for nx,ny in [(512,256),(2048,1024)]:
    for ref,test,tol in [('cpu-float','gpu-float','5e-4'),('cpu-double','gpu-double','1e-9'),('cpu-double','gpu-float','5e-4')]:
        a=f'{ref}-{nx}x{ny}.bin';b=f'{test}-{nx}x{ny}.bin'
        r=subprocess.run([str(root/'build/hos_compare'),str(out/a),str(out/b),tol],capture_output=True,text=True)
        record={'reference':a,'candidate':b,'tolerance':tol,'exit':r.returncode,
                'fields':[json.loads(line) for line in r.stdout.splitlines() if line.startswith('{')],'stderr':r.stderr}
        checks.append(record);print(json.dumps(record),flush=True)
        (out/'comparisons.json').write_text(json.dumps(checks,indent=2))
        if r.returncode:raise RuntimeError('Numerical comparison failed')
nvml.nvmlShutdown()
