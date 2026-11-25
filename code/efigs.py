import json, math, re, numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
S=json.load(open('stats.json')); R=S['res']
plt.rcParams.update({'font.family':'Liberation Sans','font.size':7.5})
COL={'Congenital/perinatal viral (CMV)':'#c0392b','Invasive CNS':'#2a78d6','Invasive systemic':'#eb6834','Respiratory':'#1baf7a','Mucocutaneous/enteric/exanthem':'#7a5fd0','Noninfectious/inflammatory':'#8a8a85','Reference':'#444'}
short={'Urinary tract infection/pyelonephritis':'UTI/pyelonephritis','Respiratory syncytial virus infection':'RSV infection','Herpes simplex virus infection':'HSV infection','Diaper dermatitis (reference)':'Diaper dermatitis'}
nm=lambda k: short.get(R[k]['name'],R[k]['name'])
ks=[k for k in R if k!='W23']
leg=[Line2D([],[],marker='s',ls='',color=c,ms=5,label=g) for g,c in COL.items() if g!='Reference']
def save(fig,name):
    for ext in ['png','pdf']: fig.savefig(f'figs/{name}.{ext}',dpi=300,bbox_inches='tight')
# eFigure: calibrated HR
fig,axs=plt.subplots(1,2,figsize=(7.2,5.2),dpi=300,sharey=False)
for ax,o,t in zip(axs,['Autism','ADHD'],['A  Autism, calibrated','B  ADHD, calibrated']):
    kk=sorted(ks,key=lambda k:R[k]['out'][o]['cal'])
    for y,k in enumerate(kk):
        v=R[k]['out'][o]; c=COL[R[k]['group']]
        ax.plot([v['cal_lo'],v['cal_hi']],[y,y],color=c,lw=1.1); ax.plot(v['cal'],y,'s',color=c,ms=4)
    ax.axvline(1,color='k',lw=0.8); ax.set_xscale('log'); ax.set_xticks([0.5,0.75,1,1.5,2]); ax.set_xticklabels(['0.5','0.75','1','1.5','2']); ax.minorticks_off()
    ax.set_yticks(range(len(kk))); ax.set_yticklabels([nm(k) for k in kk],fontsize=7)
    for lab in ax.get_yticklabels():
        if lab.get_text()=='CMV disease': lab.set_fontweight('bold')
    ax.set_xlabel('Ratio of HRs vs diaper dermatitis (95% CI)'); ax.set_title(t,loc='left',fontweight='bold',fontsize=8.5)
    for s in ['top','right']: ax.spines[s].set_visible(False)
fig.legend(handles=leg,loc='lower center',ncol=3,frameon=False,fontsize=6.8,bbox_to_anchor=(0.5,-0.06))
fig.tight_layout(); save(fig,'eFigure1_calibrated_HR')
# eFigure: risk difference per 1000
fig,axs=plt.subplots(1,2,figsize=(7.2,5.2),dpi=300)
for ax,o,t in zip(axs,['Autism','ADHD'],['A  Autism','B  ADHD']):
    kk=sorted(ks+['W23'],key=lambda k:R[k]['out'][o]['rd'])
    for y,k in enumerate(kk):
        v=R[k]['out'][o]; c=COL[R[k]['group']]
        ax.barh(y,v['rd'],color=c,height=0.7)
        ax.text(v['rd']+(0.6 if v['rd']>=0 else -0.6),y,f"{v['rd']:+.1f}",va='center',ha='left' if v['rd']>=0 else 'right',fontsize=6.3)
    ax.axvline(0,color='k',lw=0.8); ax.set_yticks(range(len(kk))); ax.set_yticklabels([nm(k) for k in kk],fontsize=7)
    lo,hi=ax.get_xlim(); ax.set_xlim(lo-4,hi+4)
    ax.set_xlabel('Risk difference per 1,000 children\n(exposed minus matched controls)'); ax.set_title(t,loc='left',fontweight='bold',fontsize=8.5)
    for s in ['top','right']: ax.spines[s].set_visible(False)
fig.legend(handles=leg+[Line2D([],[],marker='s',ls='',color='#444',ms=5,label='Reference (diaper dermatitis)')],loc='lower center',ncol=4,frameon=False,fontsize=6.5,bbox_to_anchor=(0.5,-0.07))
fig.tight_layout(); save(fig,'eFigure2_risk_difference')
# eFigure: unrestricted vs restricted
X={};cur=None
for l in open('../infx/results.txt'):
    m=re.match(r'^X(\d\d) \S+ pairs (\d+)',l)
    if m: cur='W'+m.group(1); X[cur]={'pairs':int(m.group(2))}; continue
    m=re.match(r'^(\w+) \| .*\| ([\d.]+) \(([\d.]+)-([\d.]+)\)',l)
    if m and cur: X[cur][m.group(1)]=tuple(map(float,m.group(2,3,4)))
X['W00']={'pairs':5147,'Autism':(1.45,1.21,1.73),'ADHD':(0.87,0.72,1.05)}
json.dump(X,open('unrestricted.json','w'),indent=1)
fig,axs=plt.subplots(1,2,figsize=(7.2,3.7),dpi=300)
for ax,o,t in zip(axs,['Autism','ADHD'],['A  Autism','B  ADHD']):
    for k in X:
        if k not in R or o not in X[k]: continue
        a=X[k][o]; b=R[k]['out'][o]; c=COL[R[k]['group']]
        ax.plot([a[1],a[2]],[b['hr']]*2,color=c,lw=0.6,alpha=0.6); ax.plot([a[0]]*2,[b['lo'],b['hi']],color=c,lw=0.6,alpha=0.6)
        ax.plot(a[0],b['hr'],'o',color=c,ms=4.5,mec='white',mew=0.4)
        if k in ('W00','W12','W23','W01','W03'): ax.annotate(nm(k),(a[0],b['hr']),xytext=(4,3),textcoords='offset points',fontsize=6.3)
    ax.set_xscale('log'); ax.set_yscale('log'); lim=(0.5,2.2); ax.set_xlim(lim); ax.set_ylim(lim)
    ax.plot(lim,lim,color='#999',ls='--',lw=0.7); ax.axhline(1,color='k',lw=0.5); ax.axvline(1,color='k',lw=0.5)
    tk=[0.5,0.75,1,1.5,2]; ax.set_xticks(tk); ax.set_yticks(tk); ax.set_xticklabels([f'{t:g}' for t in tk]); ax.set_yticklabels([f'{t:g}' for t in tk]); ax.minorticks_off()
    ax.set_xlabel('HR, exposure not required to have a well-child visit'); ax.set_ylabel('HR, well-child visit required (primary)'); ax.set_title(t,loc='left',fontweight='bold',fontsize=8.5)
    for s in ['top','right']: ax.spines[s].set_visible(False)
fig.legend(handles=leg+[Line2D([],[],marker='s',ls='',color='#444',ms=5,label='Reference (diaper dermatitis)')],loc='lower center',ncol=4,frameon=False,fontsize=6.5,bbox_to_anchor=(0.5,-0.1))
fig.tight_layout(); save(fig,'eFigure3_restricted_vs_unrestricted')
# eFigure: E-values for associations with CI excluding 1 (HR>1) for autism/ADHD/secondary
items=[]
for k in R:
    if k=='W23': continue
    for o in ['Autism','ADHD','Speech','ID','CP','Epilepsy','SNHL']:
        v=R[k]['out'].get(o)
        if v and v['lo']>1 and v.get('q') is not None and v['q']<0.05: items.append((nm(k),o,v['hr'],v['e_pt'],v['e_ci'],R[k]['group']))
items.sort(key=lambda x:x[3])
OL={'Autism':'autism','ADHD':'ADHD','Speech':'speech-language','ID':'intellectual disability','CP':'cerebral palsy','Epilepsy':'epilepsy','SNHL':'hearing loss'}
fig,ax=plt.subplots(figsize=(7.2,0.16*len(items)+1.2),dpi=300)
for y,it in enumerate(items):
    c=COL[it[5]]; ax.barh(y,it[3],color=c,alpha=0.35,height=0.72); ax.barh(y,it[4],color=c,height=0.72)
    ax.text(it[3]+0.08,y,f'{it[3]:.2f} / {it[4]:.2f}',va='center',fontsize=6)
ax.set_yticks(range(len(items))); ax.set_yticklabels([f'{a} → {OL[b]} (HR {h:.2f})' for a,b,h,*_ in items],fontsize=6.4)
ax.axvline(1,color='k',lw=0.6); ax.set_xlabel('E-value (light bar: point estimate; dark bar: CI limit closest to 1)')
for s in ['top','right']: ax.spines[s].set_visible(False)
ax.set_ylim(-0.6,len(items)-0.4)
fig.tight_layout(); save(fig,'eFigure4_E_values')
print(len(items))
