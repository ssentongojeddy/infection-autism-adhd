import json, math, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
S=json.load(open('stats.json')); R=S['res']; P=S['pooled']
plt.rcParams.update({'font.family':'Liberation Sans','font.size':8})
cats=[('CMV disease (single exposure)','#c0392b',R['W00']['out']),('Invasive CNS infection (k=2)','#2a78d6',P['Invasive CNS']),('Invasive systemic infection (k=3)','#eb6834',P['Invasive systemic']),('Respiratory infection (k=5)','#1baf7a',P['Respiratory']),('Mucocutaneous, enteric, exanthematous (k=6)','#7a5fd0',P['Mucocutaneous/enteric/exanthem']),('Noninfectious or inflammatory (k=4)','#8a8a85',P['Noninfectious/inflammatory'])]
outs=[('Autism','Autism'),('ADHD','ADHD'),('Speech','Speech-language\ndisorder'),('ID','Intellectual\ndisability'),('CP','Cerebral\npalsy'),('Epilepsy','Epilepsy'),('SNHL','Sensorineural\nhearing loss'),('NegCtl','Arm fracture\n(negative control)')]
fig,ax=plt.subplots(figsize=(7.2,4.4),dpi=300)
w=0.13; off=np.linspace(-2.5*w,2.5*w,len(cats))
for j,(o,_) in enumerate(outs):
    if j%2==0: ax.axvspan(j-0.5,j+0.5,color='#f3f3f1',lw=0,zorder=0)
    for (lab,c,d),dx in zip(cats,off):
        if o not in d: continue
        v=d[o]; x=j+dx
        ax.plot([x,x],[v['lo'],v['hi']],color=c,lw=1.3,zorder=2)
        sig=v['lo']>1
        ax.plot(x,v['hr'],'o',ms=5 if sig else 4.2,mfc=c if sig else 'white',mec=c,mew=1.2,zorder=3)
ax.axhline(1,color='black',lw=0.8,zorder=1)
ax.set_yscale('log'); ax.set_ylim(0.3,9)
yt=[0.33,0.5,0.75,1,1.5,2,3,5,8]; ax.set_yticks(yt); ax.set_yticklabels([f'{t:g}' for t in yt]); ax.minorticks_off()
ax.set_ylabel('Hazard ratio (95% CI), log scale')
ax.set_xticks(range(len(outs))); ax.set_xticklabels([l for _,l in outs],fontsize=7.3)
ax.set_xlim(-0.5,len(outs)-0.5)
ax.axvline(1.5,color='#999',lw=0.6,ls=':'); ax.axvline(6.5,color='#999',lw=0.6,ls=':')
ax.text(0.5,8.3,'Primary outcomes',ha='center',fontsize=7.5,fontweight='bold')
ax.text(4,8.3,'Secondary neurodevelopmental and neurologic outcomes',ha='center',fontsize=7.5,fontweight='bold')
for s in ['top','right']: ax.spines[s].set_visible(False)
h=[Line2D([],[],marker='o',ls='-',color=c,mfc=c,ms=5,lw=1.2,label=l) for l,c,_ in cats]
h+= [Line2D([],[],marker='o',ls='',mfc='white',mec='#444',ms=4.5,label='Open marker: 95% CI includes or falls below 1')]
fig.legend(handles=h,loc='lower center',ncol=3,frameon=False,fontsize=6.8,bbox_to_anchor=(0.53,0.0),columnspacing=1.2,handlelength=1.6)
fig.subplots_adjust(left=0.09,right=0.99,top=0.95,bottom=0.22)
for ext in ['png','pdf']: fig.savefig(f'figs/Figure4_risk_stratification.{ext}',dpi=300)
fig.savefig('figs/Figure4_risk_stratification.tiff',dpi=300,pil_kwargs={'compression':'tiff_lzw'})
