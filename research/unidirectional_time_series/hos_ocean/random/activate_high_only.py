import hashlib,json,pathlib,re,subprocess,time
r=pathlib.Path(__file__).resolve().parent
assert not (r/'cancel-low.json').exists()
assert 'USER_CANCELLED_LOW_BEFORE_INITIALIZATION' in (r/'prepare_random_run.m').read_text()
assert not (r/'akp002').exists(), 'Low already started: reassess before changing queue'
s=json.loads((r/'pipeline-status.json').read_text())
assert s['stage']=='HOS_running' and s['active_family']=='akp012'
config=pathlib.Path('/usr/local/bin/run.sh').read_text()
assert re.search(r'^EMAIL="\$\{EMAIL:-[^}]+\}"$',config,re.M), 'Existing mail recipient not found'
assert pathlib.Path('/usr/bin/mail').is_file()
subprocess.run(['sh','-n',str(r/'completion_email.sh')],check=True)
for name in ['finish_high_only.py','activate_high_only.py']:
    compile((r/name).read_text(),name,'exec')
utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
record={'requested_utc':utc,'action':'Cancel low before initialization; retain high HOS and GL; email success or failure once',
        'files_sha256':{name:hashlib.sha256((r/name).read_bytes()).hexdigest() for name in ['prepare_random_run.m','finish_high_only.py','completion_email.sh','activate_high_only.py']}}
(r/'cancel-low.json').write_text(json.dumps(record,indent=2))
s.update(next_family=None,low_status='cancelled_before_initialization',email='queued_after_high_HOS_and_GL')
(r/'pipeline-status.json').write_text(json.dumps(s,indent=2))
log=(r/'high-only-watcher.log').open('a')
p=subprocess.Popen(['python3',str(r/'finish_high_only.py')],cwd=r,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
(r/'high-only-watcher.pid').write_text(str(p.pid)+'\n')
print(json.dumps({'low':'cancelled_before_initialization','watcher_pid':p.pid,'email':'configured_existing_recipient_not_sent','high':'left_running'},indent=2))
