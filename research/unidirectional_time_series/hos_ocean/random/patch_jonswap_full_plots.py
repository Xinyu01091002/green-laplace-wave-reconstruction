import pathlib,hashlib,json
r=pathlib.Path(__file__).resolve().parent;p=r/'compare_jonswap_time.m';t=p.read_text();old=hashlib.sha256(p.read_bytes()).hexdigest()
a='grid on;xlim(d.limits);bound=';b="grid on;xlim([t(1) t(end)]);xline(d.limits(1),':');xline(d.limits(2),':');bound="
assert t.count(a)==1;t=t.replace(a,b);p.write_text(t)
(r/'postprocess-display-manifest.json').write_text(json.dumps({'before_sha256':old,'after_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'scope':'Display full random record; fixed interior scoring unchanged'},indent=2))