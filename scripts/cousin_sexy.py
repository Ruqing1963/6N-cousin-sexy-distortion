#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Reconnaissance (S10): omega-dependence of twin (diff 2), cousin (diff 4), and
sexy (diff 6) prime pairs on the 6N skeleton, binned by the LEFT center's
omega_{>3}. Geometry:
  twin   : (6N-1, 6N+1)            single center N, both wings
  cousin : (6N+1, 6(N+1)-1)        straddles centers N, N+1 (right wing + left wing)
  sexyA  : (6N-1, 6(N+1)-1)        straddles N, N+1, both left wings
  sexyB  : (6N+1, 6(N+1)+1)        straddles N, N+1, both right wings
Attribution: to the LEFT center's omega (same convention as the twin series).
Memory-light: streaming per-segment accumulation (does NOT store all centers).
Carries one-center overlap between segments so straddling pairs at the boundary
are counted. Default S10. Requires: numpy.
"""
import numpy as np, math, os
def primes_upto(n):
    s=np.ones(n+1,bool); s[:2]=False
    for i in range(2,int(math.isqrt(n))+1):
        if s[i]: s[i*i::i]=False
    return np.nonzero(s)[0].astype(np.int64)
MAXK=int(os.environ.get("MAXK",10))
LO=10**(MAXK-1)//6+1; HI=10**MAXK//6; SEG=4_000_000
PB=int(math.isqrt(6*HI+250))+1; BP=primes_upto(PB)
OMAX=7
# accumulators per left-center omega: total centers, and pair counts
tot=np.zeros(OMAX+2, dtype=np.int64)      # centers with given omega (for twin & as left center)
c_twin=np.zeros(OMAX+2, dtype=np.int64)
c_cous=np.zeros(OMAX+2, dtype=np.int64)
c_sexA=np.zeros(OMAX+2, dtype=np.int64)
c_sexB=np.zeros(OMAX+2, dtype=np.int64)
import time; t0=time.time()
prev_om=None; prev_pm=None; prev_pp=None; prev_lastN=None
n=LO
while n<=HI:
    nh=min(n+SEG,HI+1); sz=nh-n
    rem=np.arange(n,nh,dtype=np.int64); ob=np.zeros(sz,np.int16)
    for p in BP:
        if p*p>nh-1: break
        f=((n+p-1)//p)*p
        if f>=nh: continue
        idx=np.arange(f-n,sz,p)
        if idx.size==0: continue
        sub=rem[idx]; m=(sub%p)==0
        while m.any(): sub[m]//=p; m=(sub%p)==0
        rem[idx]=sub
        if p>3: ob[idx]+=1
    ob[rem>1]+=1
    Narr=np.arange(n,nh,dtype=np.int64)
    vlo=6*n-1; vhi=6*(nh-1)+1; span=vhi-vlo+1
    comp=np.zeros(span,bool); sq=int(math.isqrt(vhi))+1
    for p in BP:
        if p>sq: break
        st=max(p*p,((vlo+p-1)//p)*p)
        if st>vhi: continue
        comp[st-vlo:span:p]=True
    pm=~comp[(6*Narr-1)-vlo]  # 6N-1 prime
    pp=~comp[(6*Narr+1)-vlo]  # 6N+1 prime
    omc=np.clip(ob,0,OMAX+1)
    # twin: single center
    tw=pm&pp
    for om in range(1,OMAX+1):
        s=(omc==om)
        tot[om]+=s.sum()
        c_twin[om]+=(tw&s).sum()
    # straddling pairs within this segment (left center i, right center i+1)
    if sz>=2:
        leftom=omc[:-1]
        cous=pp[:-1]&pm[1:]
        sexA=pm[:-1]&pm[1:]
        sexB=pp[:-1]&pp[1:]
        for om in range(1,OMAX+1):
            s=(leftom==om)
            c_cous[om]+=(cous&s).sum()
            c_sexA[om]+=(sexA&s).sum()
            c_sexB[om]+=(sexB&s).sum()
    # boundary: pair between previous segment's last center and this segment's first
    if prev_lastN is not None and prev_lastN+1==n:
        lo_om=prev_om
        if 1<=lo_om<=OMAX:
            if prev_pp and pm[0]: c_cous[lo_om]+=1
            if prev_pm and pm[0]: c_sexA[lo_om]+=1
            if prev_pp and pp[0]: c_sexB[lo_om]+=1
    prev_om=int(omc[-1]); prev_pm=bool(pm[-1]); prev_pp=bool(pp[-1]); prev_lastN=int(Narr[-1])
    n=nh
print(f"S{MAXK}: scan {time.time()-t0:.0f}s")
print(f"\n{'omega':>5}{'Ncenters':>13}{'twin':>10}{'cousin':>10}{'sexyA':>10}{'sexyB':>10}")
for om in range(1,OMAX+1):
    if tot[om]<20000: continue
    print(f"{om:>5}{tot[om]:>13,}{c_twin[om]/tot[om]:>10.5f}{c_cous[om]/tot[om]:>10.5f}{c_sexA[om]/tot[om]:>10.5f}{c_sexB[om]/tot[om]:>10.5f}")
print(f"\nNormalised to omega=1 (shape):")
base={}
for om in range(1,OMAX+1):
    if tot[om]>=20000:
        base[om]=(c_twin[om]/tot[om],c_cous[om]/tot[om],c_sexA[om]/tot[om],c_sexB[om]/tot[om])
o1=base[1]
print(f"{'omega':>5}{'twin':>9}{'cousin':>9}{'sexyA':>9}{'sexyB':>9}")
for om in sorted(base):
    r=base[om]
    print(f"{om:>5}{r[0]/o1[0]:>9.3f}{r[1]/o1[1]:>9.3f}{r[2]/o1[2]:>9.3f}{r[3]/o1[3]:>9.3f}")

# ---- emit CSV (cousin_sexy_S{K}_data.csv) ----
import csv as _csv
with open(f'cousin_sexy_S{MAXK}_data.csv','w',newline='') as _f:
    _w=_csv.writer(_f)
    _w.writerow(['omega','Ncenters','twin','cousin','sexyA','sexyB',
                 'twin_norm','cousin_norm','sexyA_norm','sexyB_norm'])
    _b={}
    for om in range(1,OMAX+1):
        if tot[om]>=20000:
            _b[om]=(c_twin[om]/tot[om],c_cous[om]/tot[om],c_sexA[om]/tot[om],c_sexB[om]/tot[om])
    _o1=_b[min(_b)]
    for om in sorted(_b):
        r=_b[om]
        _w.writerow([om,int(tot[om]),f'{r[0]:.5f}',f'{r[1]:.5f}',f'{r[2]:.5f}',f'{r[3]:.5f}',
                     f'{r[0]/_o1[0]:.3f}',f'{r[1]/_o1[1]:.3f}',f'{r[2]/_o1[2]:.3f}',f'{r[3]/_o1[3]:.3f}'])
print(f"\n[ok] wrote cousin_sexy_S{MAXK}_data.csv")
