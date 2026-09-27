"""Temporary LAN transfer: two named files, one client IP, random capability."""
import http.server,json,pathlib,secrets,threading,sys
root=pathlib.Path(sys.argv[1]).resolve();ready=root/'reference-transfer-ready.json'
kind=sys.argv[2] if len(sys.argv)>2 else 'rhs';assert kind in ['rhs','adaptive','benchmark','campaign']
assert not ready.exists()
token=secrets.token_hex(24)
filename='adaptive_case5.bin' if kind=='adaptive' else 'ocean_case5.bin'
files={family:root/('jonswap_'+family+'_phi000/'+filename) for family in ['low','high']}
if kind=='benchmark':files={'runtime':root/'cpu-runtime.tar.gz','reference':root/'reference-records.tar.gz'}
if kind=='campaign':files={'initials':root/'campaign-initials.tar.gz'}
assert all(p.is_file() for p in files.values())
done=set();lock=threading.Lock()
class Handler(http.server.BaseHTTPRequestHandler):
    def do_GET(self):
        family=self.path.removeprefix('/'+token+'/')
        if self.client_address[0]!='192.168.2.92' or self.path!='/'+token+'/'+family or family not in files:
            self.send_error(403);return
        path=files[family];self.send_response(200);self.send_header('Content-Length',str(path.stat().st_size));self.end_headers()
        with path.open('rb') as f:
            while chunk:=f.read(1024*1024):self.wfile.write(chunk)
        self.wfile.flush()
        with lock:
            done.add(family)
            if len(done)==len(files):threading.Thread(target=server.shutdown,daemon=True).start()
    def log_message(self,*args):pass
server=http.server.ThreadingHTTPServer(('192.168.2.93',0),Handler)
info={'port':server.server_port,'token':token}
ready.write_text(json.dumps(info))
print(json.dumps(info),flush=True)
timer=threading.Timer(180,server.shutdown);timer.daemon=True;timer.start()
try:server.serve_forever()
finally:timer.cancel();server.server_close();ready.unlink(missing_ok=True)
