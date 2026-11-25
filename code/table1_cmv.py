import json, math
# label, code, in_psm, before (n1,%1,n2,%2), after (n1,%1,n2,%2)
N=(1347,3090897,1345,1345)
rows=[
('Male','M',1,747,55.457,1581797,51.176,745,55.39,765,56.877),
('Female','F',0,597,44.321,1500492,48.546,597,44.387,579,43.048),
('Hispanic or Latino','2135-2',1,367,27.246,620299,20.069,366,27.212,351,26.097),
('Not Hispanic or Latino','2186-5',1,786,58.352,1748790,56.579,785,58.364,805,59.851),
('Unknown ethnicity','UN',0,194,14.402,721808,23.353,194,14.424,189,14.052),
('White','2106-3',1,564,41.871,1554000,50.277,563,41.859,569,42.305),
('Black or African American','2054-5',1,343,25.464,565742,18.303,343,25.502,354,26.32),
('Asian','2028-9',1,102,7.572,140522,4.546,102,7.584,103,7.658),
('Native Hawaiian or Other Pacific Islander','2076-8',0,31,2.301,24433,0.79,31,2.305,10,0.743),
('American Indian or Alaska Native','1002-5',0,11,0.817,18197,0.589,11,0.818,10,0.743),
('Other race','2131-1',0,112,8.315,243510,7.878,112,8.327,138,10.26),
('Unknown race','UNK',0,184,13.66,544493,17.616,183,13.606,164,12.193),
('Short gestation/low birth weight','P07',1,364,27.023,107707,3.485,362,26.914,437,32.491),
('Slow fetal growth/fetal malnutrition','P05',1,115,8.537,14743,0.477,113,8.401,99,7.361),
('Newborn affected by maternal factors or delivery complications','P00-P04',1,80,5.939,56212,1.819,80,5.948,102,7.584),
('Birth trauma','P10-P15',1,10,0.742,9473,0.306,10,0.743,15,1.115),
('Respiratory distress of newborn','P22',1,106,7.869,15036,0.486,106,7.881,105,7.807),
('Transitory endocrine/metabolic disorders of newborn','P70-P74',1,108,8.018,12074,0.391,108,8.03,101,7.509),
('Other disturbances of cerebral status of newborn','P91',1,36,2.673,6218,0.201,36,2.677,33,2.454),
('Feeding problems of newborn','P92',0,163,12.101,63232,2.046,163,12.119,176,13.086),
('Congenital malformations of nervous system','Q00-Q07',1,173,12.843,11303,0.366,171,12.714,155,11.524),
('Congenital malformations of circulatory system','Q20-Q28',1,310,23.014,64870,2.099,309,22.974,343,25.502),
('Cleft lip and palate','Q35-Q37',0,10,0.742,4561,0.148,10,0.743,15,1.115),
('Congenital musculoskeletal malformations','Q65-Q79',0,135,10.022,105946,3.428,134,9.963,146,10.855),
('Chromosomal abnormalities','Q90-Q99',1,32,2.376,9512,0.308,31,2.305,36,2.677),
('Disorders of the immune mechanism','D80-D89',1,88,6.533,1958,0.063,86,6.394,49,3.643),
('Lack of expected normal physiological development','R62',1,301,22.346,69705,2.255,300,22.305,308,22.9),
('Symptoms concerning food and fluid intake','R63',1,304,22.569,93265,3.017,303,22.528,300,22.305),
('Acute upper respiratory infection','J06',1,197,14.625,199340,6.449,195,14.498,202,15.019),
('Asthma','J45',1,24,1.782,16641,0.538,24,1.784,25,1.859),
('Nonsuppurative otitis media','H65',1,53,3.935,28354,0.917,53,3.941,63,4.684),
('Suppurative/unspecified otitis media','H66',1,96,7.127,88933,2.877,95,7.063,100,7.435),
('Gastroesophageal reflux disease','K21',0,190,14.105,156747,5.071,188,13.978,221,16.431),
('Family history of mental and behavioral disorders','Z81',1,10,0.742,3633,0.118,10,0.743,10,0.743),
('Family history of disabilities and chronic diseases','Z82',0,45,3.341,17567,0.568,45,3.346,35,2.602),
('Well-child visit (any time to index)','Z00.1',0,756,56.125,3090897,100.0,755,56.134,1345,100.0),
('Systemic antibacterials','ATC J01',1,441,32.739,231558,7.492,440,32.714,499,37.1),
]
def smd(p1,p2):
    p1/=100;p2/=100; d=math.sqrt((p1*(1-p1)+p2*(1-p2))/2); return abs(p1-p2)/d if d>0 else 0
def smdc(m1,s1,m2,s2): return abs(m1-m2)/math.sqrt((s1**2+s2**2)/2)
out=[]
out.append(dict(label='Current age, mean (SD), y',code='',psm=1,b1='6.66 (4.47)',b2='7.45 (4.81)',a1='6.66 (4.48)',a2='6.65 (4.54)',sb=smdc(6.66,4.47,7.45,4.81),sa=smdc(6.66,4.48,6.65,4.54)))
out.append(dict(label='Age at index, mean (SD), y',code='',psm=1,b1='0.18 (0.38)',b2='0.16 (0.36)',a1='0.18 (0.38)',a2='0.18 (0.39)',sb=smdc(.179,.383,.155,.362),sa=smdc(.178,.383,.184,.387)))
for r in rows:
    l,c,p,n1,p1,n2,p2,a1,q1,a2,q2=r
    f=lambda n,pc: ('≤10' if n<=10 else f'{n:,}')+f' ({pc:.1f})'
    out.append(dict(label=l,code=c,psm=p,b1=f(n1,p1),b2=f(n2,p2),a1=f(a1,q1),a2=f(a2,q2),sb=smd(p1,p2),sa=smd(q1,q2)))
json.dump(dict(N=N,rows=out),open('table1_cmv.json','w'),indent=1)
for o in out: print(f"{o['label'][:40]:40} {o['b1']:>14} {o['b2']:>18} {o['sb']:.3f} | {o['a1']:>12} {o['a2']:>12} {o['sa']:.3f} {'PSM' if o['psm'] else ''}")
