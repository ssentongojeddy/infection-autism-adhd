import json, math
S=json.load(open('stats.json')); R=S['res']; P=S['pooled']; X=json.load(open('unrestricted.json')); T=json.load(open('table1_cmv.json'))
order=['W00','W01','W02','W06','W09','W08','W10','W11','W12','W13','W14','W03','W04','W15','W16','W17','W18','W07','W19','W21','W22','W23']
short={'Urinary tract infection/pyelonephritis':'UTI/pyelonephritis','Respiratory syncytial virus infection':'RSV infection','Herpes simplex virus infection':'HSV infection'}
nm=lambda k: short.get(R[k]['name'],R[k]['name'])
OUT=[('Autism','Autism spectrum disorder'),('ADHD','ADHD'),('Speech','Speech-language disorder'),('ID','Intellectual disability'),('CP','Cerebral palsy'),('Epilepsy','Epilepsy'),('SNHL','Sensorineural hearing loss'),('NegCtl','Arm fracture (negative-control outcome)')]
wq={'W01':1348,'W02':1859,'W03':13795,'W04':40698,'W06':10654,'W07':1881,'W08':60499,'W09':1381,'W10':166721,'W11':142086,'W12':95708,'W13':119589,'W14':3354,'W15':124821,'W16':12410,'W17':5954,'W18':154545,'W19':32660,'W21':396999,'W22':66097,'W23':502903,'W00':1387}
cq={'W%02d'%i:None for i in range(24)}
for l in open('../infx/counts.txt'):
    p=l.split()
    if p[0].startswith('C'): cq['W'+p[0][1:]]=int(p[1])
cq['W00']=3228691
xc={}
for l in open('../infx/counts.txt'):
    p=l.split()
    if p[0].startswith('X'): xc['W'+p[0][1:]]=(p[2],int(p[3]))
def fmt(n): return f'{n:,}'
def f2(x): return f'{x:.2f}'
def pf(p):
    if p<0.001: return '<.001'
    return f'{p:.3f}'.lstrip('0') if p<0.01 else f'{p:.2f}'.lstrip('0')
def ev(n): return '≤10' if n<=10 else fmt(n)
md={}
# Table 1 main (condensed)
rows=T['rows']; lines=['| Characteristic | CMV disease, before (n = 1,347) | Controls, before (n = 3,090,897) | SMD, before | CMV disease, after (n = 1,345) | Controls, after (n = 1,345) | SMD, after |','|:--|--:|--:|--:|--:|--:|--:|']
for r in rows:
    lab=r['label']+(f" ({r['code']})" if r['code'] and not r['code'][0].isdigit() and r['code'] not in('M','F','UN','UNK') else '')
    if not r['psm']: lab+=' ^a^'
    lines.append(f"| {lab} | {r['b1']} | {r['b2']} | {r['sb']:.2f} | {r['a1']} | {r['a2']} | {r['sa']:.2f} |")
md['TABLE1']='\n'.join(lines)+'\n\nAbbreviations: CMV, cytomegalovirus; SMD, absolute standardized mean difference.\n\nValues are No. (%) unless otherwise indicated. Before-matching denominators are infants with data in the baseline window as reported by TriNetX (cohort sizes 1,386 and 3,228,691). TriNetX reports counts of 10 or fewer as 10; these are shown as ≤10, and SMDs for them are approximate. Race and ethnicity were as recorded in the electronic health record.\n\n^a^ Not included in the propensity score model (shown to assess balance on variables not used for matching). The well-child visit was required for all infants before age 2 years; in the exposed group it could occur after the CMV index date, so it appears in the covariate window for only 56.1%.'
json.dump(md,open('main_tables.json','w'))
# Supplement tables
sup={}
# eTable1 exposures
L=['| Group | Condition | ICD-10-CM codes | Infants with code before age 1 y, No. ^a^ | With well-child visit, No. ^b^ | Eligible controls, No. | Matched pairs, No. |','|:--|:--|:--|--:|--:|--:|--:|']
G={'Congenital/perinatal viral (CMV)':'Congenital/perinatal viral','Reference':'Benign reference'}
for k in order:
    a=xc.get(k,('',None))[1]; 
    if k=='W00': a=5518
    L.append(f"| {G.get(R[k]['group'],R[k]['group'])} | {nm(k)} | {R[k]['code']} | {fmt(a) if a else '—'} | {fmt(wq[k])} | {fmt(cq[k]) if cq.get(k) else '—'} | {fmt(R[k]['pairs'])} |")
sup['eT1']='\n'.join(L)+'\n\n^a^ First code before age 2 years, without the well-child requirement (from the cohort query at the time of analysis). ^b^ Primary exposure cohort: code and a well-child visit (Z00.1) before age 2 years. Planned but not analyzable because TriNetX restricts the codes: congenital CMV (P35.1), neonatal sepsis (P36), and neonatal jaundice (P59).'
# eTable 4 per outcome full results
blocks=[]
for o,ol in OUT:
    L=[f'**{ol}**','','| Condition | At risk, exposed | Events, exposed | At risk, controls | Events, controls | Risk/1000, exposed | Risk/1000, controls | RD/1000 | HR (95% CI) | *P* | FDR *q* | E-value (CI) |','|:--|--:|--:|--:|--:|--:|--:|--:|:--|--:|--:|--:|']
    for k in order:
        v=R[k]['out'].get(o)
        if not v: L.append(f'| {nm(k)} | NR | NR | NR | NR | | | | NR | | | |'); continue
        q='—' if v.get('q') is None else pf(v['q'])
        L.append(f"| {nm(k)} | {fmt(v['n1'])} | {ev(v['e1'])} | {fmt(v['n2'])} | {ev(v['e2'])} | {v['risk1']:.1f} | {v['risk2']:.1f} | {v['rd']:+.1f} | {f2(v['hr'])} ({f2(v['lo'])}-{f2(v['hi'])}) | {pf(v['p'])} | {q} | {v['e_pt']:.2f} ({v['e_ci']:.2f}) |")
    blocks.append('\n'.join(L))
sup['eT4']='\n\n'.join(blocks)
# eTable 5 calibrated
L=['| Condition | '+' | '.join(l for o,l in OUT[:-1])+' |','|:--|'+'--:|'*7]
L[0]='| Condition | Autism | ADHD | Speech-language | Intellectual disability | Cerebral palsy | Epilepsy | Hearing loss |'
for k in order[:-1]:
    c=[]
    for o,_ in OUT[:-1]:
        v=R[k]['out'].get(o)
        c.append(f"{f2(v['cal'])} ({f2(v['cal_lo'])}-{f2(v['cal_hi'])})" if v and 'cal' in v else 'NR')
    L.append(f"| {nm(k)} | "+' | '.join(c)+' |')
ref=R['W23']['out']
L.append('| *Reference HR, diaper dermatitis* | '+' | '.join(f"*{f2(ref[o]['hr'])} ({f2(ref[o]['lo'])}-{f2(ref[o]['hi'])})*" for o,_ in OUT[:-1])+' |')
sup['eT5']='\n'.join(L)
# eTable 6 pooled
L=['| Condition group (k) | Outcome | Pooled HR (95% CI) | *I*^2^, % |','|:--|:--|:--|--:|']
kk={'Invasive CNS':2,'Invasive systemic':3,'Respiratory':5,'Mucocutaneous/enteric/exanthem':6,'Noninfectious/inflammatory':4}
for g,d in P.items():
    for o,ol in OUT:
        v=d.get(o)
        if v: L.append(f"| {g} ({kk[g]}) | {ol} | {f2(v['hr'])} ({f2(v['lo'])}-{f2(v['hi'])}) | {v['I2']:.0f} |")
sup['eT6']='\n'.join(L)
# eTable 7 unrestricted (long format)
def h(t): return f'{f2(t[0])} ({f2(t[1])}-{f2(t[2])})'
def hv(v): return f"{f2(v['hr'])} ({f2(v['lo'])}-{f2(v['hi'])})"
L=['| Condition | Matched pairs, unrestricted / primary | Outcome | HR (95% CI), unrestricted | HR (95% CI), primary |','|:--|:--|:--|:--|:--|']
ON={'Autism':'Autism','ADHD':'ADHD','Speech':'Speech-language disorder','ID':'Intellectual disability','CP':'Cerebral palsy','Epilepsy':'Epilepsy','SNHL':'Hearing loss','NegCtl':'Arm fracture'}
for k in order:
    if k not in X: continue
    x=X[k]; r=R[k]['out']; first=True
    for o in ON:
        if o not in x: continue
        L.append(f"| {nm(k) if first else ''} | {fmt(x['pairs'])+' / '+fmt(R[k]['pairs']) if first else ''} | {ON[o]} | {h(x[o])} | {hv(r[o]) if o in r else 'NR'} |"); first=False
sup['eT7']='\n'.join(L)+'\n\nUnrestricted: exposure cohort did not require a well-child visit (controls did). Bacterial meningitis, viral meningitis/encephalitis, and HSV infection used 30 covariates (adding P92, K21, and Z82); others used the same 27 covariates as the primary analysis. CMV unrestricted estimates are from a separate earlier analysis matched on 9 covariates (age at index, current age, sex, race, ethnicity, P07, Q20-Q28, D80-D89, Q90-Q99). Unrestricted analyses were not completed for febrile seizure, atopic dermatitis, and iron-deficiency anemia.'
# eTable 8 E-values for significant positive associations
items=[]
for k in order[:-1]:
    for o,ol in OUT[:-1]:
        v=R[k]['out'].get(o)
        if v and v['lo']>1 and v.get('q') is not None and v['q']<0.05: items.append((v['e_pt'],nm(k),ol,v))
items.sort(reverse=True)
L=['| Condition | Outcome | HR (95% CI) | FDR *q* | E-value, point | E-value, CI limit |','|:--|:--|:--|--:|--:|--:|']
for e,n,ol,v in items: L.append(f"| {n} | {ol} | {hv(v)} | {pf(v['q'])} | {v['e_pt']:.2f} | {v['e_ci']:.2f} |")
sup['eT8']='\n'.join(L)
# eTable 9 full CMV baseline (same as table1 but complete with codes)
sup['eT9']=md['TABLE1']
json.dump(sup,open('supp_tables.json','w'))
print('ok',len(items))
