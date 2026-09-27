"""On 93: export four frozen low-family initial states through official modules."""
import hashlib,json,pathlib,subprocess,sys,tarfile
r=pathlib.Path(__file__).resolve().parent
p=r/'ocean_adaptive_probe.f90';s=p.read_text();assert s.count('T_stop_star=.3_rp/T')==1
p.write_text(s.replace('T_stop_star=.3_rp/T','T_stop_star=0.0_rp'))
s=(r/'build_ocean_probe.py').read_text()
s=s.replace("for family in ['low','high']:","for family,phase in [('low',p) for p in [0,90,180,270]]:")
s=s.replace("origin=campaign/family/'cases/phi000';case=root/('jonswap_'+family+'_phi000')","origin=campaign/family/('cases/phi%03d'%phase);case=root/('jonswap_'+family+('_phi%03d'%phase))")
s=s.replace("'phase':0,'source_parts'","'phase':phase,'source_parts'")
assert "origin=campaign/family/'cases/phi000'" not in s
(r/'campaign_builder.py').write_text(s)
subprocess.run([sys.executable,str(r/'campaign_builder.py'),'--production','--adaptive'],check=True)
origin=pathlib.Path('/home/lxy/green-laplace-unidirectional-time-series-runs/hos-jonswap-r4gl-80tp-20260926-v2/low')
files=[]
for phase in [0,90,180,270]:
    name='phi%03d'%phase;initial=r/('jonswap_low_'+name)/'adaptive_case5.bin'
    target=r/(name+'.bin');target.write_bytes(initial.read_bytes());assert target.stat().st_size==108+2*513*512*16
    lines=(origin/'cases'/name/'Results/probes.dat').read_text().splitlines();i=next(i for i,s in enumerate(lines) if 'VARIABLES' in s)
    row=next(s.split() for s in lines[i+1:] if s.strip());assert float(row[0])==0 and len(row)==6
    ref=r/(name+'-initial.csv');ref.write_text('time,p1,p2,p3,p4,p5,potential,kinetic,total,relative_energy_change\n'+','.join(row+['0']*4)+'\n')
    files.extend([target,ref])
files.extend([r/'actual_input_provenance.json',r/'provenance.json'])
(r/'settings.json').write_bytes((origin/'settings.json').read_bytes());files.append(r/'settings.json')
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in files};(r/'initial-manifest.json').write_text(json.dumps(manifest,indent=2));files.append(r/'initial-manifest.json')
with tarfile.open(r/'campaign-initials.tar.gz','w:gz') as tar:
    for p in files:tar.add(p,arcname=p.name)
print(json.dumps({'bundle_sha256':hashlib.sha256((r/'campaign-initials.tar.gz').read_bytes()).hexdigest()}),flush=True)
