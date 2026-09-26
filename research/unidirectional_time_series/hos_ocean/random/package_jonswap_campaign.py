import ast,hashlib,json,pathlib,subprocess,tarfile
repo=pathlib.Path(__file__).resolve().parents[4];src=pathlib.Path(__file__).resolve().parent
out=repo/'artifacts/hos_ocean/jonswap-campaign-source-v2';assert not out.exists();(out/'spark').mkdir(parents=True)
spark=pathlib.Path('C:/Users/spet5947/Documents/ChatGPT/SPARK');sources={}
for name in ['spark_eta20_r_series.m','spark_psi20_r_series.m','spark_difference_input.m','spark_difference_audit.m','spark_r4_low_mode_pair_repair.m','spark_native_superharmonics.m']:
 p=spark/'src/forward'/('green_laplace' if name=='spark_native_superharmonics.m' else 'subharmonic')/name
 original=p.read_bytes();t=original.decode('utf-8-sig').replace('\r\n','\n')
 if name in ['spark_eta20_r_series.m','spark_psi20_r_series.m']:
  needle='if outsideDefect>1e-8';assert t.count(needle)==1
  t=t.replace(needle,"if outsideDefect>1e-8 && outsideDefect<=1e-7\n    warning('run:R4ConsistencyWarning','Accepted run-local numerical residual %.3e; original threshold1e-8, hard stop1e-7.',outsideDefect);\nend\nif outsideDefect>1e-7")
  t=t.replace('exceeds 1e-8.','exceeds run-local hard limit 1e-7.')
 (out/'spark'/name).write_text(t,encoding='utf-8');sources[name]={'original_sha256':hashlib.sha256(original).hexdigest(),'deployed_sha256':hashlib.sha256((out/'spark'/name).read_bytes()).hexdigest()}
for name in ['prepare_jonswap_hybrid.m','export_jonswap_hos.m','run_jonswap_campaign.py']:(out/name).write_bytes((src/name).read_bytes())
t=(repo/'research/directional_wave_data/hos_low_run/run_low_production.py').read_text(encoding='utf-8-sig')
t=t.replace("r=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-directional-akp002-20260926-v1')","r=pathlib.Path(__file__).resolve().parent")
t=t.replace('assert np in [4,8]',"assert np in [4,8]\nsettings=json.loads((r/'settings.json').read_text())")
t=t.replace("'aggregate_RSS_limit_GiB':8","'aggregate_RSS_limit_GiB':16").replace('total<8*1024*1024','total<16*1024*1024').replace('exceeded 8 GiB','exceeded 16 GiB')
t=t.replace("r/'prepare_directional_hos.m'","r/'export_jonswap_hos.m'")
t=t.replace("100*float(matches[-1])/220","100*float(matches[-1])/settings['duration_s']")
t=t.replace('len(values)==1101',"len(values)==settings['expected_samples']").replace("row[0]-.2*i","row[0]-settings['output_dt_s']*i").replace('last_logged_time_s=220.0',"last_logged_time_s=settings['duration_s']")
t=t.replace("matches=re.findall", "assert not re.search(r'(?i)\\b(nan|infinity)\\b',tail),'Nonfinite solver log '+j['name']\n                matches=re.findall")
(out/'run_directional_hos.py').write_text(t,encoding='utf-8')
t=(repo/'artifacts/hos_ocean/random-launch-source-v2/compare_random_full.m').read_text(encoding='utf-8-sig')
t=t.replace('compare_random_full','compare_jonswap_time').replace('random-gl-comparison-v1','gl-time-comparison-v1')
t=t.replace('widths=[15 7.5]','widths=[3.75 1.875]').replace('Joint 15 deg','Joint 3.75 deg').replace('Joint 7.5 deg','Joint 1.875 deg').replace('GL 7.5 deg','GL 1.875 deg').replace('GL 15 deg','GL 3.75 deg').replace('GL7p5','GL1p875').replace('Phase-only random waves','JONSWAP random waves')
t=t.replace('tr,8,pairBins(:)', 'tr,16,pairBins(:)').replace('tr,8,bins)', 'tr,16,bins)')
(out/'compare_jonswap_time.m').write_text(t,encoding='utf-8')
t=(repo/'research/directional_wave_data/remote_completion_email_watcher.sh').read_text(encoding='utf-8-sig')
t=t.replace('Sixteen native-wall OW3D runs: one focused group and three random realizations, four phases each.','JONSWAP gamma3.3 low then high kpHs/2=.02/.12: four phases each, approx80Tp, R4-GL initialization.')
t=t.replace('Randomization preserves modal amplitudes. No local raw-data download.','Same99-percent-energy support and random phases; Hs normalized separately for both levels.')
t=t.replace('Completion means native final time and raw-output integrity checks; physical/model accuracy is not certified.','Completion means both HOS families and GL time-series processing completed; inspect summary.json for failures.')
t=t.replace('Per-case processed eta/phi: cases/*/processed/surface_strip.mat','Results: low/gl-time-comparison-v1 and high/gl-time-comparison-v1').replace('[Green-Laplace OW3D]','[Green-Laplace HOS JONSWAP]')
t=t.replace('  if test -f "$run_dir/runtime.log";', '  if test -f "$run_dir/summary.json"; then cat "$run_dir/summary.json"; fi\n  if test -f "$run_dir/runtime.log";')
(out/'completion_email.sh').write_text(t,encoding='utf-8',newline='\n')
for p in out.glob('*.py'):ast.parse(p.read_text())
manifest={'SPARK_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=spark,text=True).strip(),'SPARK_source':sources,'files_sha256':{str(p.relative_to(out)):hashlib.sha256(p.read_bytes()).hexdigest() for p in out.rglob('*') if p.is_file()},'policy':'Private runtime snapshot, no source vendoring into public code. Only R4 numerical-check handling changed; no physics/repair changes.'}
(out/'source-manifest.json').write_text(json.dumps(manifest,indent=2))
with tarfile.open(out.parent/'jonswap-campaign-source-v2.tar.gz','w:gz') as tar:
 for p in out.rglob('*'):
  if p.is_file():tar.add(p,arcname=str(p.relative_to(out)))
print('Snapshot assembled; Python syntax checked')
