#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Build the 2-panel cousin/sexy figure from ../data/cousin_sexy_S10_data.csv
(produced by cousin_sexy.py with default MAXK=10). Left: the three omega-distortions
(twin rises, cousin/sexy-A fall and overlap, sexy-B non-monotone). Right: cousin
and sexy-A curves coincide (shared right member 6(N+1)-1).
"""
import csv, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows=list(csv.DictReader(open('../data/cousin_sexy_S10_data.csv')))
om=np.array([int(r['omega']) for r in rows])
tw=np.array([float(r['twin_norm']) for r in rows])
co=np.array([float(r['cousin_norm']) for r in rows])
sa=np.array([float(r['sexyA_norm']) for r in rows])
sb=np.array([float(r['sexyB_norm']) for r in rows])
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13.5,5.4))
ax1.plot(om,tw,'o-',color='#c0392b',lw=2.4,ms=8,label='twin (diff 2, single-centre)',zorder=5)
ax1.plot(om,co,'s-',color='#185FA5',lw=2,ms=7,label='cousin (diff 4)',zorder=4)
ax1.plot(om,sa,'^--',color='#3aa0d0',lw=1.6,ms=7,label='sexy-A (diff 6, left wings)',zorder=3)
ax1.plot(om,sb,'D-',color='#2ca25f',lw=2,ms=7,label='sexy-B (diff 6, right wings)',zorder=4)
ax1.axhline(1,color='gray',ls=':',lw=1)
ax1.set_xlabel(r'$\omega_{>3}$ of left centre $N$',fontsize=11)
ax1.set_ylabel(r'pair rate, normalised to $\omega=1$',fontsize=11)
ax1.set_title('Three distinct $\\omega$-distortions on the $6N$ skeleton',fontsize=12)
ax1.legend(fontsize=9,loc='upper left'); ax1.grid(alpha=.25); ax1.set_xticks(range(1,8))
ax2.plot(om,co,'s-',color='#185FA5',lw=2,ms=9,label='cousin (diff 4)',zorder=3)
ax2.plot(om,sa,'^',color='#e8845b',ms=9,mew=2,label='sexy-A (diff 6, left wings)',zorder=4)
ax2.set_xlabel(r'$\omega_{>3}$ of left centre $N$',fontsize=11)
ax2.set_ylabel(r'pair rate, normalised',fontsize=11)
ax2.set_title('cousin $\\equiv$ sexy-A: identical $\\omega$-curves\n(shared right member $6(N{+}1){-}1$)',fontsize=12)
ax2.legend(fontsize=10); ax2.grid(alpha=.25); ax2.set_xticks(range(1,8))
plt.suptitle('Cousin and sexy prime pairs in $S_{10}$ (15B centres): $\\omega$-dependence distorts by pairing geometry',fontsize=12.5,y=1.02)
plt.tight_layout()
plt.savefig('fig_paper9_cousin_sexy.pdf',bbox_inches='tight')
plt.savefig('fig_paper9_cousin_sexy.png',dpi=160,bbox_inches='tight')
print("figure saved")
