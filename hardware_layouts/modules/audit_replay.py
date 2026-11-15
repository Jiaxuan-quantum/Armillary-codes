"""Reconstruct every physical pairing directly from the frozen code parameters.

No geometry generator is used to obtain the required edge union. This is a
circuit and code identity audit, not a new exhaustive distance calculation.
"""
from pathlib import Path
import json, hashlib, sys
from collections import Counter
from route_model import read,write

HERE=Path(__file__).resolve().parent

def replay(record,model):
    q=model['q'];p=record['p'];M=12;P=p*q;n=M*P
    channel=record['mode']=='channel'
    def site(r,c,i,j):return (r*M+c)*P+q*i+j
    def unpack(s):
        g,m=divmod(s,P);r,c=divmod(g,M);i,j=divmod(m,q)
        return r,c,i,j
    at=[0]*(2*n);pos=[0]*(2*n)
    for a in range(M):
        for i in range(p):
            for j in range(q):
                anc=a*P+q*i+j;dat=n+anc
                for lab,r in ((anc,1),(dat,0)):
                    s=site(r,(-a)%M,i,j);at[s]=lab;pos[lab]=s
    edge_union=set();counts=Counter();steps=[]
    def permutation(t,a):
        top=record['top']
        if isinstance(top,dict):return top['permutations'][f'{t},{a}']
        f,r=top
        return [((-1)**f[t*M+a]*(j-r[t*M+a]))%q for j in range(q)]
    xrows=set(record['skeleton']['xrows']);zrows=set(record['skeleton']['zrows'])
    orbit_keys=set();step_pairs=[]
    for t in range(M):
        pairs=[];covered=set()
        for a in range(M):
            perm=permutation(t,a) if t%2==0 else None
            for m in range(P):
                s=pos[a*P+m];r,c,i,j=unpack(s)
                if t%2==0:
                    assert r==1;k=perm[j];v=site(0,c,i,k)
                    key=('v',c,k,j)
                else:
                    assert r==0;d=record['phases'][(t//2*M+a)*(q if channel else 1)+(j if channel else 0)]
                    v=site(1,(c-1)%M,(i+d)%p,j);key=('h',c,j,d)
                assert s not in covered and v not in covered and at[v]>=n
                covered.update((s,v));pairs.append((s,v));edge_union.add(tuple(sorted((s,v))))
                orbit_keys.add(key)
                counts['CXSWAP_X' if a in xrows else 'CXSWAP_Z' if a in zrows else 'SWAP_only']+=1
        assert len(covered)==2*n
        for s,v in pairs:
            aa,bb=at[s],at[v];at[s],at[v]=bb,aa;pos[aa],pos[bb]=v,s
        steps.append(dict(step=t,operations=len(pairs),every_physical_position_used_once=True))
        step_pairs.append(pairs)
    exported=set()
    for orbit in model['orbits']:
        r,c,iz,j=orbit['source'];s,cc,jz,k=orbit['target']
        for ii in range(p):
            exported.add(tuple(sorted((site(r,c,(iz+ii)%p,j),site(s,cc,(jz+ii)%p,k)))))
    assert edge_union==exported
    assert orbit_keys=={tuple(e['key']) for e in model['orbits']}
    assert len(edge_union)==model['total_connections'] and len(orbit_keys)==len(model['orbits'])
    radii=sorted(d['radius'] for d in model['devices'])
    assert len(set(radii))==24*q
    import math
    a0=min(min(b-a for a,b in zip(radii,radii[1:])),2*min(radii)*math.sin(math.pi/p))
    assert abs(a0-1.)<1e-10
    digest=lambda x:hashlib.sha256(json.dumps(x,separators=(',',':')).encode()).hexdigest()
    return dict(status='PASS',code=model['code'],total_positions=2*n,steps=steps,
                operations=dict(counts),total_operations=sum(counts.values()),
                edge_orbits=len(orbit_keys),unique_physical_edges=len(edge_union),
                exported_full_edge_union_matches_replayed_schedule=True,
                edge_union_sha256=digest(sorted(edge_union)),schedule_sha256=digest(step_pairs),
                minimum_device_center_distance=a0,
                circuit_change='None. Geometry and face assignments are outside the fixed circuit pairing and gate schedule.',
                evidence_scope='Distance inherited from frozen code identity and archived evidence. Neither a new distance exhaustion nor a noise simulation is performed.')
