"""Sequential longer GPU runs, bounded and monitored without device-setting writes."""
import ctypes,json,pathlib,subprocess,time
root=pathlib.Path(__file__).resolve().parent
initial=pathlib.Path('/root/hos-production-check-20260927T025257Z')
lib=ctypes.CDLL('libnvidia-ml.so.1');assert lib.nvmlInit_v2()==0
device=ctypes.c_void_p();assert lib.nvmlDeviceGetHandleByIndex_v2(0,ctypes.byref(device))==0
class Memory(ctypes.Structure):_fields_=[('total',ctypes.c_ulonglong),('free',ctypes.c_ulonglong),('used',ctypes.c_ulonglong)]
def memory():
    x=Memory();assert lib.nvmlDeviceGetMemoryInfo(device,ctypes.byref(x))==0;return x.used
def temperature():
    x=ctypes.c_uint();assert lib.nvmlDeviceGetTemperature(device,0,ctypes.byref(x))==0;return x.value
rows=[]
for family in ['low','high']:
    for precision,tol in [('double','1e-12'),('float','1e-9')]:
        name=family+'-'+precision;out=root/name;assert not out.exists()
        cmd=[str(root/'build/hos_run'),precision,str(initial/(family+'-adaptive.bin')),str(out),'27.6',tol]
        baseline=memory();peak=baseline;temp=temperature();start=time.monotonic();last=start
        print('Starting '+name+' for 27.6 physical seconds',flush=True)
        with (root/(name+'.log')).open('w') as log:
            process=subprocess.Popen(cmd,stdout=log,stderr=subprocess.STDOUT)
            while process.poll() is None:
                peak=max(peak,memory());temp=max(temp,temperature())
                if time.monotonic()-last>=30:
                    try:status=json.loads((out/'status.json').read_text())
                    except (OSError,json.JSONDecodeError):status={}
                    print(json.dumps({'run':name,'progress':status,'elapsed_seconds':time.monotonic()-start}),flush=True);last=time.monotonic()
                if time.monotonic()-start>3600:
                    process.terminate();process.wait();raise TimeoutError(name)
                time.sleep(.1)
        row={'run':name,'exit':process.returncode,'gpu_baseline_bytes':baseline,'gpu_peak_bytes':peak,'gpu_increment_peak_bytes':peak-baseline,'temperature_peak_sampled':temp,'memory_sampling_seconds':.1}
        if (out/'summary.json').exists():row['summary']=json.loads((out/'summary.json').read_text())
        else:row['error']=(root/(name+'.log')).read_text()[-2000:]
        rows.append(row);(root/'matrix.json').write_text(json.dumps(rows,indent=2));print(json.dumps(row),flush=True)
lib.nvmlShutdown()
