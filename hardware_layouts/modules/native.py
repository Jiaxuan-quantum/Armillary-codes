"""Independent physical-position reconstruction for the frozen q3/q4 rules.

No old generator, inverse synthesis, geometry optimizer, or motif producer is
imported. All removed checks still take part in the twelve transport steps.
"""
from collections import Counter, defaultdict
import hashlib, json


def perm_at(record, q, t, a):
    top = record['top']
    if isinstance(top, dict):
        return top['permutations'][f'{t},{a}']
    f, r = top
    return [((-1) ** f[t * 12 + a] * (j-r[t * 12 + a])) % q for j in range(q)]


def physical_forms(record, q, channel):
    nv = 72 * (q if channel else 1)
    at = {}; pos = {}; coords = {}
    for a in range(12):
        for j in range(q):
            for rail in (0, 1):
                label = (rail, a, j); site = (rail, (-a) % 12, j)
                at[site] = label; pos[label] = site; coords[label] = [0] * nv
    blocks = {}
    for t in range(12):
        pairs = []; covered = set()
        for a in range(12):
            for qa in range(q):
                label = (1, a, qa); src = pos[label]; rail, c, j = src
                delta = [0] * nv
                if t % 2 == 0:
                    assert rail == 1
                    dst = (0, c, perm_at(record, q, t, a)[j])
                else:
                    assert rail == 0
                    dst = (1, (c-1) % 12, j)
                    delta[(t//2*12+a)*(q if channel else 1)+(j if channel else 0)] = 1
                assert src not in covered and dst not in covered
                covered.update((src, dst)); other = at[dst]
                assert other[0] == 0 and other[1] == (a+t) % 12
                blocks[a, other[1], qa] = (other[2], [x+d-y for x,d,y in zip(coords[label],delta,coords[other])])
                pairs.append((src, dst, label, other, delta))
        assert len(covered) == 24*q
        for src, dst, label, other, delta in pairs:
            at[src], at[dst] = other, label; pos[label], pos[other] = dst, src
            coords[label] = [x+d for x,d in zip(coords[label],delta)]
            coords[other] = [x-d for x,d in zip(coords[other],delta)]
    return blocks


def echelon(rows):
    basis = {}
    for r in rows:
        while r:
            pivot = r.bit_length()-1
            if pivot in basis:
                r ^= basis[pivot]
            else:
                basis[pivot] = r; break
    return basis


def remainder(v, basis):
    while v:
        p = v.bit_length()-1
        if p not in basis:
            break
        v ^= basis[p]
    return v


def verify_logical(code, witness):
    side = witness['basis'].lower(); support = witness['support']
    assert side in ('x','z') and len(set(support)) == len(support)
    assert all(0 <= k < code['n'] for k in support)
    v = sum(1 << k for k in support)
    assert all((v&r).bit_count()%2 == 0 for r in code['hz' if side=='x' else 'hx'])
    assert remainder(v, echelon(code['h'+side])) != 0
    return len(support)


def literal_code(record, q, phases=None, channel=None):
    if phases is None:
        phases = record['phases']
    if channel is None:
        channel = record['mode'] == 'channel'
    p = record['p']; P = q*p; n = 12*P
    def site(rail,c,z,j):
        return (rail*12+c)*P+q*z+j
    at = [0]*(2*n); pos = [0]*(2*n)
    for a in range(12):
        for z in range(p):
            for j in range(q):
                anc = a*P+q*z+j; dat = n+anc
                s = site(1,(-a)%12,z,j); at[s]=anc; pos[anc]=s
                s = site(0,(-a)%12,z,j); at[s]=dat; pos[dat]=s
    initial = at[:]; steps=[]
    for t in range(12):
        pairs=[]; covered=set()
        for a in range(12):
            perm = perm_at(record,q,t,a) if t%2==0 else None
            for member in range(P):
                label=a*P+member; src=pos[label]
                group,w=divmod(src,P);rail,c=divmod(group,12);z,j=divmod(w,q)
                if t%2==0:
                    assert rail==1; dst=site(0,c,z,perm[j])
                else:
                    assert rail==0
                    delta=phases[(t//2*12+a)*(q if channel else 1)+(j if channel else 0)]
                    dst=site(1,(c-1)%12,(z+delta)%p,j)
                assert src not in covered and dst not in covered and at[dst]>=n
                covered.update((src,dst));pairs.append((src,dst))
        assert len(covered)==2*n
        for u,v in pairs:
            a,b=at[u],at[v];at[u],at[v]=b,a;pos[a],pos[b]=v,u
        steps.append(pairs)
    at=initial[:];both=[];masks=[record['skeleton']['xrows'],record['skeleton']['zrows']]
    xs,zs=map(set,masks)
    for sequence in (steps,list(reversed(steps))):
        bits=[1<<a for a in at];vis=defaultdict(list)
        for pairs in sequence:
            for u,v in pairs:
                aa,dd=(u,v) if at[u]<n else (v,u)
                a,d=at[aa],at[dd]-n
                assert 0<=a<n and 0<=d<n
                vis[a].append(d)
                if a//P in xs:bits[dd]^=bits[aa]
                elif a//P in zs:bits[aa]^=bits[dd]
                at[u],at[v]=at[v],at[u];bits[u],bits[v]=bits[v],bits[u]
        incidence=defaultdict(int)
        for a,ds in vis.items():
            assert len(ds)==len(set(ds))==12 and sorted(d//P for d in ds)==list(range(12))
            if a//P in xs:
                for d in ds:incidence[d]^=1<<a
        for s,label in enumerate(at):
            expected=1<<label
            if label>=n:expected^=incidence[label-n]
            elif label//P in zs:expected^=sum(1<<(n+d) for d in vis[label])
            assert bits[s]==expected
        both.append([[sorted(vis[a*P+j]) for a in mask for j in range(P)] for mask in masks])
    assert at==initial and both[0]==both[1]
    hx,hz=[[sum(1<<d for d in row) for row in rows] for rows in both[0]]
    girth=True
    for rows in both[0]:
        pairs=Counter((row[i],row[j]) for row in rows for i in range(12) for j in range(i+1,12))
        girth = girth and max(pairs.values())==1
    digest=hashlib.sha256(json.dumps([n,list(map(hex,hx)),list(map(hex,hz))],separators=(',',':')).encode()).hexdigest()
    rx,rz=len(echelon(hx)),len(echelon(hz))
    return dict(n=n,k=n-rx-rz,hx=hx,hz=hz,rank_X=rx,rank_Z=rz,matrix_sha256=digest,
                girth_ge6=girth,native_both_rounds_literal=True,slots_restored=True)
