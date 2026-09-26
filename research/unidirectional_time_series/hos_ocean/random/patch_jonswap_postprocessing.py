"""Drop exactly zero directional coefficients before unchanged pair evaluation."""
import hashlib,json,pathlib
r=pathlib.Path(__file__).resolve().parent
p=r/'compare_jonswap_time.m';t=p.read_text();old=hashlib.sha256(p.read_bytes()).hexdigest()
assert 'function [e,p]=nonzero_gl_sum' not in t
t=t.replace('gl_directional_sum_time(','nonzero_gl_sum(')
t+='''
function [e,p]=nonzero_gl_sum(A,w,kx,ky,g,h,kp,t,J,bins)
active=A(:)~=0;
assert(any(active));
[e,p]=gl_directional_sum_time(A(active),w(active),kx(active),ky(active),g,h,kp,t,J,bins(active));
end
'''
p.write_text(t)
(r/'postprocess-runtime-patch.json').write_text(json.dumps({'before_sha256':old,'after_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'change':'Remove exactly zero amplitudes only, no tolerance cutoff or kernel change','email_script_sha256':hashlib.sha256((r/'completion_email.sh').read_bytes()).hexdigest()},indent=2))
