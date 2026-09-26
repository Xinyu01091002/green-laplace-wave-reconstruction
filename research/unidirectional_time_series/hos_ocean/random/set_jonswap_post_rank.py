import pathlib,hashlib,json
r=pathlib.Path(__file__).resolve().parent;p=r/'compare_jonswap_time.m';t=p.read_text();old=hashlib.sha256(p.read_bytes()).hexdigest()
assert 'tr,8,pairBins(:)' in t
t=t.replace('tr,8,pairBins(:)','tr,16,pairBins(:)').replace('tr,8,bins)','tr,16,bins)')
p.write_text(t)
(r/'postprocess-rank-manifest.json').write_text(json.dumps({'before_sha256':old,'after_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'eta22_GL_rank':16,'reason':'Use same selected GL rank for broad-band eta22 time reconstruction as initialization; no kernel change'},indent=2))