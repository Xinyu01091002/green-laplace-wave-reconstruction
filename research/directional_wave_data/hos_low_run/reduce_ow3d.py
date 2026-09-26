import pathlib,json,os,signal,time
r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/ow3d16-kpd1-akp012-20260925T212822Z');s=json.loads((r/'status.json').read_text());controller=int((r/'controller.pid').read_text());targets=[];kept=[]
for name,v in s['cases'].items():
 if name.startswith(('random_s20260926_','random_s20260927_')):targets.append((name,v))
 else:kept.append((name,v))
assert len(targets)==8 and len(kept)==8
for name,v in targets+kept:
 p=pathlib.Path('/proc',str(v['pid']));assert (p/'comm').read_text().strip()=='ow3d';assert (p/'cwd').resolve()==(r/'cases'/name).resolve();assert (p/'exe').resolve()==(r/'bin/ow3d').resolve();assert os.getpgid(v['pid'])==v['pid'];assert int((p/'stat').read_text().split(') ')[1].split()[1])==controller
note={'requested_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),'reason':'User requested retain one focused wavegroup plus one random realization, each four phases; release CPU/RAM for HOS','retained':[x[0] for x in kept],'stopped':[{'case':n,'pid':v['pid']} for n,v in targets],'files_preserved':True,'controller_not_stopped':True,'note':'Legacy frozen controller may label intentionally terminated cases failed; use this user-requested cancellation record to interpret them.'}
p=r/'user-requested-reduction-20260926.json';assert not p.exists();p.write_text(json.dumps(note,indent=2))
for name,v in targets:
 (r/'cases'/name/'user-cancelled.json').write_text(json.dumps({'reason':note['reason'],'utc':note['requested_utc'],'signal':'SIGTERM'},indent=2));os.killpg(v['pid'],signal.SIGTERM)
time.sleep(3)
note['retained_live']=[n for n,v in kept if pathlib.Path('/proc',str(v['pid'])).exists()];note['stopped_still_live']=[n for n,v in targets if pathlib.Path('/proc',str(v['pid'])).exists()];assert len(note['retained_live'])==8;assert not note['stopped_still_live'];p.write_text(json.dumps(note,indent=2));print(json.dumps(note,indent=2))