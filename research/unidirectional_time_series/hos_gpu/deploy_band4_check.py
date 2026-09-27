"""Replace the narrow edge-band pilot with the user-requested 0<omega<=4omega_p test."""
import hashlib,json,pathlib,shlex,subprocess
s=pathlib.Path(__file__).resolve().parent;repo=s.parents[2];private=repo/'artifacts/hos_gpu/higher-harmonics-v1';local=repo/'artifacts/hos_gpu/higher-20260927T145332Z/first-band-4wp-v1';local.mkdir(exist_ok=False)
r='/home/lxy/green-laplace-unidirectional-time-series-runs/gpu-higher-harmonics-20260927T145332Z'
body=(private/'check_first_harmonic_band.m').read_text()
body=body.replace('check_first_harmonic_band(run)','check_first_harmonic_band4(run)').replace('first-band-check-v1','first-band-4wp-v1')
old="variants={baseBins,(min(baseBins)-2:max(baseBins)+2)',(lowAllowed:ceil(1.1*max(baseBins)))'};names={'original_band','two_bins_each_side','wider_admissible_band'};"
assert body.count(old)==1;body=body.replace(old,"variants={baseBins,allBins(allW<=4*2*pi/s.Tp)};names={'original_band','all_nonzero_to_4wp'};")
body=body.replace('zeros(N,3)','zeros(N,numel(variants))').replace('cell(1,3)','cell(1,numel(variants))').replace('for j=1:3','for j=1:numel(variants)')
body=body.replace('assert(all(kk*s.h>.5)&&3*max(ww)<pi/dt)','assert(all(kk>0)&&3*max(ww)<pi/dt)')
body=body.replace('gl_third_directional_time(','gl_third_directional_time_fullband(').replace("legend(['Unprojected',names]","legend([{'Unprojected'},names]")
oldnote='Controlled finite-record spectral-support sensitivity, not a claim that out-of-initial-band FFT leakage is physical free-wave energy. Original in-band directional coefficients unchanged. Raw, fixed-common-band and fixed-taper reference definitions remain fixed across variants. No truly unfiltered all-frequency GL run: current third-order domain and directional information do not justify treating all observed FFT bins as free parents.'
newnote='User-requested input test: retain every nonzero temporal FFT bin at or below 4omega_p, without original-spectrum min/max or kh>.5 input removal. Original in-band coefficients and all output/reference processing remain fixed. Low-kh quadrature convergence and the finite-record directional prior outside its physical source band are unverified; this is a diagnostic, not physical certification of all admitted components. No gain or time alignment.'
assert oldnote in body;body=body.replace(oldnote,newnote);m=local/'check_first_harmonic_band4.m';m.write_text(body)
kernel=(private/'gl_third_directional_time.m').read_text();assert kernel.count('all(q>.5)')==1
kernel=kernel.replace('gl_third_directional_time(','gl_third_directional_time_fullband(').replace('all(q>.5)','all(q>0)')
k=local/'gl_third_directional_time_fullband.m';k.write_text(kernel)
controller=(s/'run_band_controller.py').read_text().replace('first-band-check-v1','first-band-4wp-v1').replace('check_first_harmonic_band(','check_first_harmonic_band4(')
c=local/'run_band4_controller.py';c.write_text(controller)
opts=['-i',str(pathlib.Path.home()/'.ssh/id_ed25519_cursor'),'-o','BatchMode=yes'];ssh=['ssh','-p','60093',*opts,'root@60.188.112.99']
stop=r'''
import json,os,pathlib,signal,time
r=pathlib.Path(ROOT);old=r/'first-band-check-v1';s=json.loads((old/'status.json').read_text())
if s['state']=='running':
 pid=s['pid'];p=pathlib.Path('/proc',str(pid),'cmdline')
 if p.exists():
  args=p.read_bytes();assert b'check_first_harmonic_band' in args and str(r).encode() in args;assert os.getpgid(pid)==pid
  os.killpg(pid,signal.SIGTERM)
(old/'superseded.json').write_text(json.dumps({'reason':'User requested all nonzero first-harmonic bins through 4omega_p instead of narrow edge-padding trials. Existing partial results retained.','previous_status':s},indent=2))
'''.replace('ROOT',repr(r))
subprocess.run([*ssh,'python3 -c '+shlex.quote(stop)],check=True)
subprocess.run(['scp','-P','60093',*opts,str(m),str(k),str(c),'root@60.188.112.99:'+r+'/'],check=True)
hashes={f.name:hashlib.sha256(f.read_bytes()).hexdigest() for f in [m,k,c]}
code="import pathlib,subprocess; r=pathlib.Path("+repr(r)+"); assert not (r/'first-band-4wp-v1').exists(); f=(r/'band4-controller.log').open('x'); p=subprocess.Popen(['runuser','-u','lxy','--','python3',str(r/'run_band4_controller.py')],cwd=r,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True); print(p.pid)"
p=subprocess.run([*ssh,'python3 -c '+shlex.quote(code)],text=True,capture_output=True,check=True);print(p.stdout);(local/'launch.json').write_text(json.dumps({'remote':r+'/first-band-4wp-v1','pid':int(p.stdout.strip()),'source_hashes':hashes},indent=2))
