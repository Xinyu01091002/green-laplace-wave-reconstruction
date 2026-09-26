"""Reconcile the live legacy queue's intentional low cancellation and notify."""
import json
import pathlib
import subprocess
import time

root=pathlib.Path(__file__).resolve().parent
def utc(): return time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())
def write(name,value):
    p=root/name;q=p.with_suffix(p.suffix+'.tmp');q.write_text(value);q.replace(p)

write('high-only-watcher-status.txt','waiting_for_high_HOS_and_GL\n')
while True:
    state=json.loads((root/'pipeline-status.json').read_text())
    if state['stage'] in ['failed','completed_HOS_and_GL','HOS_complete_GL_needs_attention']:
        break
    time.sleep(15)
high=state.get('outcomes',{}).get('akp012',{})
success=high.get('HOS')=='completed' and high.get('GL')=='completed'
low=root/'akp002'
assert not (low/'cases').exists(), 'Unexpected low production: investigate before claiming cancellation'
summary={'state':'completed_high_only' if success else 'high_only_needs_attention',
         'updated_utc':utc(),'high':high,'low':'cancelled_by_user_before_initialization',
         'original_pipeline_terminal_state':state,'scope':'Only phase-only random Akp focusing label .12, seed20260925'}
if success:
    assert state['stage']=='failed' and 'USER_CANCELLED_LOW_BEFORE_INITIALIZATION' in state.get('error',''), 'Unexpected pipeline terminal state'
    summary['HOS_resources']=json.loads((root/'akp012/status.json').read_text())
    summary['GL_report']=str(root/'akp012/random-gl-comparison-v1/report.json')
write('summary.json',json.dumps(summary,indent=2))
write('pipeline-status.json',json.dumps(dict(stage=summary['state'],updated_utc=utc(),outcomes={'akp012':high,'akp002':{'HOS':'cancelled_by_user_not_run'}}),indent=2))
write('status.txt',summary['state']+'\n')
write('exit-code.txt',('0' if success else '1')+'\n')
write('finished.utc',utc()+'\n')
write('high-only-watcher-status.txt','finished_waiting\n')
subprocess.run(['sh',str(root/'completion_email.sh'),str(root),'hos-random-high-only-20260926'],check=True)
