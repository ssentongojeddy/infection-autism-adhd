import re, math, json
from scipy.stats import norm
names={'W00':'CMV disease','W01':'Bacterial meningitis','W02':'Viral meningitis/encephalitis','W03':'Herpes simplex virus infection','W04':'Enterovirus infection','W06':'Sepsis (postneonatal)','W07':'Kawasaki disease','W08':'Urinary tract infection/pyelonephritis','W09':'Osteomyelitis/septic arthritis','W10':'Respiratory syncytial virus infection','W11':'Influenza','W12':'COVID-19','W13':'Pneumonia','W14':'Pertussis','W15':'Viral gastroenteritis','W16':'Roseola','W17':'Varicella','W18':'Hand, foot, and mouth disease','W19':'Febrile seizure','W21':'Atopic dermatitis','W22':'Iron-deficiency anemia','W23':'Diaper dermatitis (reference)'}
codes={'W00':'B25','W01':'G00','W02':'A87, A86','W03':'B00','W04':'B34.1, B97.1','W06':'A40, A41','W07':'M30.3','W08':'N39.0, N10','W09':'M86, M00','W10':'J21.0, J12.1, B97.4','W11':'J09-J11','W12':'U07.1','W13':'J13, J15, J18','W14':'A37','W15':'A08','W16':'B08.2','W17':'B01','W18':'B08.4','W19':'R56.0','W21':'L20','W22':'D50','W23':'L22'}
grp={'W00':'Congenital/perinatal viral (CMV)','W01':'Invasive CNS','W02':'Invasive CNS','W06':'Invasive systemic','W09':'Invasive systemic','W08':'Invasive systemic','W10':'Respiratory','W11':'Respiratory','W12':'Respiratory','W13':'Respiratory','W14':'Respiratory','W03':'Mucocutaneous/enteric/exanthem','W04':'Mucocutaneous/enteric/exanthem','W15':'Mucocutaneous/enteric/exanthem','W16':'Mucocutaneous/enteric/exanthem','W17':'Mucocutaneous/enteric/exanthem','W18':'Mucocutaneous/enteric/exanthem','W07':'Noninfectious/inflammatory','W19':'Noninfectious/inflammatory','W21':'Noninfectious/inflammatory','W22':'Noninfectious/inflammatory','W23':'Reference'}
D={};cur=None
for l in open('../infx/results_w.txt'):
    m=re.match(r'^(W\d\d) \S+ pairs (\d+)',l)
    if m: cur=m.group(1); D[cur]={'pairs':int(m.group(2)),'out':{}}; continue
    m=re.match(r'^(\w+) \| ([\d,]+)/([\d,]+) vs ([\d,]+)/([\d,]+) \| ([\d.]+) \(([\d.]+)-([\d.]+)\)',l)
    if m and cur:
        n1,e1,n2,e2=[int(x.replace(',','')) for x in m.group(2,3,4,5)]
        h,lo,hi=map(float,m.group(6,7,8))
        D[cur]['out'][m.group(1)]=dict(n1=n1,e1=e1,n2=n2,e2=e2,hr=h,lo=lo,hi=hi)
def se(o): return (math.log(o['hi'])-math.log(o['lo']))/3.92
def pval(o): z=math.log(o['hr'])/se(o); return 2*(1-norm.cdf(abs(z)))
def evalue(h,lo,hi):
    def ev(r):
        if r<1: r=1/r
        return r+math.sqrt(r*(r-1))
    e_pt=ev(h)
    if lo<=1<=hi: e_ci=1.0
    else: e_ci=ev(lo if h>1 else hi)
    return e_pt,e_ci
outs=['Autism','ADHD','Speech','ID','CP','Epilepsy','SNHL','NegCtl']
# BH FDR across all exposure (excluding ref) x 7 non-negctl outcomes
tests=[]
for k in D:
    if k=='W23': continue
    for o in outs[:-1]:
        if o in D[k]['out']: tests.append((k,o,pval(D[k]['out'][o])))
tests.sort(key=lambda x:x[2]); m=len(tests); q={}
prev=1
for i in range(m-1,-1,-1):
    k,o,p=tests[i]; val=min(prev,p*m/(i+1)); prev=val; q[(k,o)]=val
ref=D['W23']['out']
res={}
for k in D:
    res[k]={'name':names[k],'code':codes[k],'group':grp[k],'pairs':D[k]['pairs'],'out':{}}
    for o,v in D[k]['out'].items():
        r1=v['e1']/v['n1']*1000; r2=v['e2']/v['n2']*1000
        e=evalue(v['hr'],v['lo'],v['hi'])
        d=dict(v); d.update(risk1=r1,risk2=r2,rd=r1-r2,p=pval(v),q=q.get((k,o)),e_pt=e[0],e_ci=e[1])
        if k!='W23' and o in ref:
            R=ref[o]; lr=math.log(v['hr'])-math.log(R['hr']); s=math.sqrt(se(v)**2+se(R)**2)
            d.update(cal=math.exp(lr),cal_lo=math.exp(lr-1.96*s),cal_hi=math.exp(lr+1.96*s))
        res[k]['out'][o]=d
# pooled random effects (DerSimonian-Laird) by group
pooled={}
for g in ['Invasive CNS','Invasive systemic','Respiratory','Mucocutaneous/enteric/exanthem','Noninfectious/inflammatory']:
    pooled[g]={}
    for o in outs:
        ys=[];ws=[]
        for k in D:
            if grp[k]==g and o in D[k]['out'] and k not in ('W00',):  # exclude CMV from CNS pool; reported separately
                v=D[k]['out'][o]; ys.append(math.log(v['hr'])); ws.append(1/se(v)**2)
        if len(ys)<2: continue
        W=sum(ws); ybar=sum(w*y for w,y in zip(ws,ys))/W; Q=sum(w*(y-ybar)**2 for w,y in zip(ws,ys)); df=len(ys)-1
        C=W-sum(w*w for w in ws)/W; t2=max(0,(Q-df)/C)
        wr=[1/(1/w+t2) for w in ws]; Wr=sum(wr); yr=sum(w*y for w,y in zip(wr,ys))/Wr; sr=math.sqrt(1/Wr)
        I2=max(0,(Q-df)/Q)*100 if Q>0 else 0
        pooled[g][o]=dict(hr=math.exp(yr),lo=math.exp(yr-1.96*sr),hi=math.exp(yr+1.96*sr),k=len(ys),I2=I2)
json.dump({'res':res,'pooled':pooled},open('stats.json','w'),indent=1)
for k in res:
    a=res[k]['out']['Autism'];b=res[k]['out']['ADHD']
    print(k,res[k]['name'][:22].ljust(22),res[k]['pairs'],'AUT %.2f q=%.3g E=%.2f/%.2f risk %.1f vs %.1f cal %s | ADHD %.2f q=%.3g cal %s'%(a['hr'],a['q'] or -1,a['e_pt'],a['e_ci'],a['risk1'],a['risk2'],('%.2f'%a['cal']) if 'cal' in a else '-',b['hr'],b['q'] or -1,('%.2f (%.2f-%.2f)'%(b['cal'],b['cal_lo'],b['cal_hi'])) if 'cal' in b else '-'))

for g in pooled:
  print(g)
  for o,v in pooled[g].items(): print('   ',o,'%.2f (%.2f-%.2f) k=%d I2=%.0f'%(v['hr'],v['lo'],v['hi'],v['k'],v['I2']))
