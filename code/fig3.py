import json, math, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
S=json.load(open('stats.json')); R=S['res']
plt.rcParams.update({'font.family':'Liberation Sans','font.size':8})
order=[('Congenital viral',['W00']),('Invasive CNS infection',['W01','W02']),('Invasive systemic infection',['W06','W09','W08']),('Respiratory infection',['W10','W11','W12','W13','W14']),('Mucocutaneous, enteric, and exanthematous infection',['W03','W04','W15','W16','W17','W18']),('Noninfectious or inflammatory condition',['W07','W19','W21','W22']),('Benign reference condition',['W23'])]
cols=['Autism','ADHD','Speech','ID','CP','Epilepsy','SNHL','NegCtl']
clab=['Autism','ADHD','Speech-language\ndisorder','Intellectual\ndisability','Cerebral palsy','Epilepsy','Sensorineural\nhearing loss','Arm fracture\n(negative control)']
short={'Urinary tract infection/pyelonephritis':'UTI/pyelonephritis','Respiratory syncytial virus infection':'RSV infection','Herpes simplex virus infection':'HSV infection','Diaper dermatitis (reference)':'Diaper dermatitis'}
rows=[]
for g,ks in order:
    rows.append(('H',g)); rows+= [('R',k) for k in ks]
n=len(rows)
cmap=LinearSegmentedColormap.from_list('d',['#1c5cab','#86b6ef','#f4f4f2','#f2a07b','#b8341e'])
norm=TwoSlopeNorm(vmin=-1.5,vcenter=0,vmax=2.5)
fig=plt.figure(figsize=(7.2,9.0),dpi=300)
ax=fig.add_axes([0.30,0.09,0.60,0.78])
ax.set_xlim(-0.5,len(cols)-0.5); ax.set_ylim(n-0.5,-0.5)
for i,(t,k) in enumerate(rows):
    if t=='H':
        ax.text(-0.5-3.6,i+0.1,k,ha='left',va='center',fontsize=7.6,fontweight='bold',color='#222',clip_on=False); continue
    for j,o in enumerate(cols):
        x=j+(0.12 if j>=2 else 0)+(0.12 if j>=7 else 0)
        if o not in R[k]['out']:
            ax.add_patch(plt.Rectangle((x-0.5,i-0.5),1,1,facecolor='#ececec',hatch='////',edgecolor='#c4c4c4',lw=0)); ax.text(x,i,'NR',ha='center',va='center',fontsize=6.5,color='#555'); continue
        v=R[k]['out'][o]; lg=math.log2(v['hr'])
        ax.add_patch(plt.Rectangle((x-0.5,i-0.5),1,1,facecolor=cmap(norm(lg)),edgecolor='white',lw=0.6))
        star='*' if (k!='W23' and v.get('q') is not None and v['q']<0.05) else ''
        ax.text(x,i,f"{v['hr']:.2f}{star}",ha='center',va='center',fontsize=7,color='white' if abs(lg)>1.1 else '#111',fontweight='bold' if (v['lo']>1 and k!='W23') else 'normal')
    ax.text(-0.5-3.45,i,short.get(R[k]['name'],R[k]['name']),ha='left',va='center',fontsize=7.5,clip_on=False,fontweight='bold' if k=='W00' else 'normal')
ax.set_xlim(-0.5,len(cols)-0.5+0.24)
xs=[j+(0.12 if j>=2 else 0)+(0.12 if j>=7 else 0) for j in range(len(cols))]
for x,l in zip(xs,clab): ax.text(x,-0.9,l,rotation=45,ha='left',va='bottom',fontsize=7.3,rotation_mode='anchor')
ax.axis('off')
cax=fig.add_axes([0.36,0.05,0.50,0.012])
g=np.linspace(-1.5,2.5,400); cax.imshow(g[None,:],cmap=cmap,norm=norm,aspect='auto',extent=(-1.5,2.5,0,1))
ticks=[0.35,0.5,0.75,1,1.5,2,3,5]; cax.set_xticks([math.log2(t) for t in ticks]); cax.set_xticklabels([f'{t:g}' for t in ticks],fontsize=7); cax.set_yticks([])
for sp in cax.spines.values(): sp.set_linewidth(0.4)
cax.set_xlabel('Hazard ratio vs propensity-matched well-child controls (log scale)',fontsize=7.5)
for ext in ['png','pdf']: fig.savefig(f'figs/Figure3_heatmap.{ext}',dpi=300)
fig.savefig('figs/Figure3_heatmap.tiff',dpi=300,pil_kwargs={'compression':'tiff_lzw'})
