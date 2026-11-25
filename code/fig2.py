import json, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, matplotlib.ticker as mt
from matplotlib.lines import Line2D
S=json.load(open('stats.json')); R=S['res']
plt.rcParams.update({'font.family':'Liberation Sans','font.size':8})
COL={'Congenital/perinatal viral (CMV)':'#c0392b','Invasive CNS':'#2a78d6','Invasive systemic':'#eb6834','Respiratory':'#1baf7a','Mucocutaneous/enteric/exanthem':'#7a5fd0','Noninfectious/inflammatory':'#8a8a85'}
short={'Urinary tract infection/pyelonephritis':'UTI/pyelonephritis','Respiratory syncytial virus infection':'RSV infection','Herpes simplex virus infection':'HSV infection'}
def panel(fig,rect,o,label,title):
    keys=[k for k in R if k!='W23' and o in R[k]['out']]
    keys.sort(key=lambda k:-R[k]['out'][o]['hr'])
    n=len(keys)
    axT=fig.add_axes(rect); axT.set_xlim(0,1); axT.set_ylim(-2.3,n+1.7); axT.axis('off')
    px0,pw=0.555,0.27
    ax=fig.add_axes([rect[0]+px0*rect[2],rect[1],pw*rect[2],rect[3]],sharey=axT)
    ax.set_facecolor('none')
    ref=R['W23']['out'][o]
    ax.fill_betweenx([-0.5,n-0.5],ref['lo'],ref['hi'],color='#cfcfcf',alpha=0.7,lw=0)
    ax.vlines(ref['hr'],-0.5,n-0.5,color='#6b6b6b',ls=(0,(2,2)),lw=0.9)
    ax.vlines(1,-0.5,n-0.5,color='black',lw=0.8)
    for i,k in enumerate(keys):
        v=R[k]['out'][o]; y=n-1-i; c=COL[R[k]['group']]
        ax.plot([v['lo'],v['hi']],[y,y],color=c,lw=1.2,solid_capstyle='butt')
        ax.plot(v['hr'],y,'s',color=c,ms=3.5+min(3.5,v['n1']/60000),mec='white',mew=0.4)
        b='bold' if k=='W00' else 'normal'
        axT.text(0.0,y,short.get(R[k]['name'],R[k]['name']),va='center',fontsize=7.4,fontweight=b)
        axT.text(0.275,y,f"{R[k]['pairs']:,}",va='center',ha='right',fontsize=7.2)
        axT.text(0.40,y,f"{v['e1']:,} ({v['risk1']/10:.1f})",va='center',ha='right',fontsize=7.2)
        axT.text(0.53,y,f"{v['e2']:,} ({v['risk2']/10:.1f})",va='center',ha='right',fontsize=7.2)
        axT.text(1.0,y,f"{v['hr']:.2f} ({v['lo']:.2f}-{v['hi']:.2f})",va='center',ha='right',fontsize=7.2,fontweight=b)
    for x,t,ha in [(0.0,'Early-life condition','left'),(0.275,'Matched\npairs, No.','right'),(0.40,'Events, exposed,\nNo. (%)','right'),(0.53,'Events, controls,\nNo. (%)','right'),(1.0,'HR (95% CI)','right')]:
        axT.text(x,n-0.15,t,fontsize=7.2,fontweight='bold',ha=ha,va='bottom')
    axT.text(0,n+1.15,f'{label}  {title}',fontsize=9.5,fontweight='bold',va='bottom')
    axT.plot([0,1],[n-0.35,n-0.35],color='black',lw=0.6)
    ax.set_xscale('log'); ax.set_xlim(0.45,2.6); ax.set_xticks([0.5,0.75,1,1.5,2.5])
    ax.xaxis.set_major_formatter(mt.FuncFormatter(lambda v,_:f'{v:g}')); ax.xaxis.set_minor_formatter(mt.NullFormatter())
    for s in ['top','right','left']: ax.spines[s].set_visible(False)
    ax.spines['bottom'].set_position(('data',-0.7)); ax.tick_params(axis='y',left=False,labelleft=False); ax.tick_params(axis='x',labelsize=7,length=2.5)
    ax.text(1.0,-2.0,'HR (95% CI), log scale',ha='center',va='top',fontsize=7.2,transform=ax.get_xaxis_transform() if False else ax.transData)
fig=plt.figure(figsize=(7.2,9.4),dpi=300)
panel(fig,[0.03,0.535,0.95,0.43],'Autism','A','Autism spectrum disorder')
panel(fig,[0.03,0.065,0.95,0.43],'ADHD','B','Attention-deficit/hyperactivity disorder (ADHD)')
fig.legend(handles=[Line2D([],[],marker='s',ls='',color=c,ms=5,label=g) for g,c in COL.items()]+[Line2D([],[],ls=(0,(2,2)),color='#6b6b6b',lw=1,label='Reference condition (diaper dermatitis): HR, with 95% CI shaded')],loc='lower center',ncol=3,frameon=False,fontsize=7,bbox_to_anchor=(0.5,0.0))
for ext in ['png','pdf']: fig.savefig(f'figs/Figure2_ranked_forest.{ext}',dpi=300)
fig.savefig('figs/Figure2_ranked_forest.tiff',dpi=300,pil_kwargs={'compression':'tiff_lzw'})
