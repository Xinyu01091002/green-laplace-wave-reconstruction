"""Plot read-only nvidia-smi telemetry; this does not model wave physics."""
import csv,datetime,json,pathlib,statistics,sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
p=pathlib.Path(sys.argv[1]);rows=[]
for row in csv.reader(p.read_text(encoding='utf-8-sig').splitlines()):
    if row and row[0].startswith('2026/'):
        rows.append([v.strip() for v in row])
times=[datetime.datetime.strptime(r[0],'%Y/%m/%d %H:%M:%S.%f') for r in rows]
x=[(t-times[0]).total_seconds() for t in times]
power=[float(r[2].split()[0]) for r in rows];temp=[float(r[1]) for r in rows]
metrics={'samples':len(rows),'start_gpu_timestamp':rows[0][0],'end_gpu_timestamp':rows[-1][0],'duration_s':x[-1],
 'power_W':{'min':min(power),'max':max(power),'mean':statistics.mean(power)},'temperature_C':{'min':min(temp),'max':max(temp)},
 'sm_clock_MHz':sorted(set(r[6] for r in rows)),'memory_clock_MHz':sorted(set(r[7] for r in rows)),'pstates':sorted(set(r[8] for r in rows)),
 'gpu_utilization':sorted(set(r[5] for r in rows))}
p.with_suffix('.json').write_text(json.dumps(metrics,indent=2))
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False})
fig,axes=plt.subplots(2,1,figsize=(10,6),sharex=True,layout='constrained')
axes[0].plot(x,power,color='#2365a6',lw=1.7,label='Measured board power')
axes[0].axhline(250,color='#9a5656',ls='--',lw=1,label='Reported power limit: 250 W')
axes[0].set(ylabel='Power (W)',ylim=(0,275));axes[0].legend(loc='center right')
axes[0].text(.02,.82,f'Mean {statistics.mean(power):.2f} W | range {min(power):.2f}–{max(power):.2f} W',transform=axes[0].transAxes)
axes[1].plot(x,temp,color='#da7827',lw=2,label='GPU temperature')
axes[1].axhline(92,color='#9a5656',ls='--',lw=1,label='Reported HW slowdown: 92 °C')
axes[1].set(ylabel='Temperature (°C)',xlabel='Elapsed time (s)',ylim=(50,100));axes[1].legend(loc='center right')
for ax in axes:ax.grid(alpha=.2);ax.set_xlim(0,max(x))
fig.suptitle('Tesla P40 — four concurrent FP64 HOS phases\n'+rows[0][0]+' to '+rows[-1][0]+' (GPU host timestamps)',fontsize=13)
fig.savefig(p.with_suffix('.png'),dpi=180);fig.savefig(p.with_suffix('.svg'))
print(json.dumps(metrics,indent=2))
