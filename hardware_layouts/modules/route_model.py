"""Explicit continuous polar routes, including all leads, for frozen layouts.

A point (u,r) means (r*cos(u/R), r*sin(u/R)). A straight segment in
these coordinates is a circular arc, radial line, or linear polar spiral.
Representative routes expand by u -> u + i*pitch, i=0..p-1.
"""
from pathlib import Path
import ast, json, math
from collections import Counter
from scipy.integrate import quad

RHO=[0,2,4,6,8,10,11,9,7,5,3,1]
ORIENTATION=[1]*6+[-1]*6

def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8-sig'))

def write(path,obj):
    Path(path).write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n',encoding='utf8')

def build(geometry,hardware):
    g=geometry; h=hardware; par=g['parameters']
    q=g['q'];p=g['p'];code=g['code']
    R=g.get('inner_radius',g.get('R')); pitch=par['pitch']
    dr=par['terminal_r_half'];du=par['terminal_u_half'];lane=par['lane_spacing']
    rho=g.get('rho',RHO);ori=g.get('orientation',ORIENTATION)
    low={};high={};rr=R
    widths={ast.literal_eval(k):v for k,v in g['local_widths'].items()}
    channel_order=g.get('physical_module_channel_order',[list(range(q)) for _ in range(12)])
    for c in sorted(range(12),key=lambda c:rho[c]):
        for j in channel_order[c]:
            low[c,j]=rr;high[c,j]=rr+widths[c,j];rr=high[c,j]+1.
    if 'low_radii' in g:
        low={ast.literal_eval(k):v for k,v in g['low_radii'].items()}
        high=({ast.literal_eval(k):v for k,v in g['high_radii'].items()}
              if 'high_radii' in g else {k:low[k]+widths[k] for k in low})
    def radius(r,c,j):
        if r==0:return high[c,j] if ori[c]>0 else low[c,j]
        c=(c+1)%12
        return low[c,j] if ori[c]>0 else high[c,j]
    groups={ast.literal_eval(k):v for k,v in g['groups'].items()}
    routes=[]
    for e in h['edge_orbits']:
        kind,c,j,k=e['key'];ra=radius(*e['source']);rb=radius(*e['target'])
        direction=1 if rb>ra else -1
        if kind=='h':
            group=next(ds for ds in groups[c,j] if k in ds)
            eta=par.get('pair_offset',.25)
            alpha=0. if len(group)==1 else (-eta if k==min(group) else eta)
            alpha=g.get('horizontal_alpha_by_key',{}).get(str(tuple(e['key'])),alpha)
            assert len(group)==1 or (len(group)==2 and group[1]-group[0]==1)
            points=[[0,ra],[alpha*pitch,ra],[alpha*pitch,ra+direction*dr],
                    [(k-alpha)*pitch,rb-direction*dr],[(k-alpha)*pitch,rb],[k*pitch,rb]]
            endpoints=[[0,ra],[k*pitch,rb]]
            source=[*e['source'][:2],0,e['source'][2]]
            target=[*e['target'][:2],k%p,e['target'][2]]
        else:
            target_ids=(sorted({x['key'][3] for x in h['edge_orbits'] if x['key'][:3]==e['key'][:3]})
                        if g.get('vertical_rank_mode','all_targets')=='actual_targets' else range(q))
            rank=sorted(target_ids,key=lambda z:abs(radius(1,c,z)-ra)).index(k)
            offset=direction*(g.get('vertical_offset_center',1.5)-rank)*lane
            sign=g.get('radial_angular_sign_by_source_channel',[(-1)**z for z in range(q)])[j]
            sign=g.get('vertical_sign_override',{}).get(str(tuple(e['key'])),sign)
            u=sign*(du+(rank+1)*lane)
            explicit=g.get('vertical_route_parameters',{}).get(str(tuple(e['key'])))
            if explicit is not None:offset=explicit['radial_offset'];u=explicit['angular_offset']
            rp=ra+offset
            points=[[0,ra],[0,rp],[u,rp],[u,rb],[0,rb]]
            endpoints=[[0,ra],[0,rb]]
            source=[*e['source'][:2],0,e['source'][2]]
            target=[*e['target'][:2],0,e['target'][2]]
        # Remove zero segments, preserving the exact remaining piecewise curve.
        clean=[points[0]]
        for point in points[1:]:
            if math.dist(point,clean[-1])>1e-13:clean.append(point)
        segment_lengths=[]
        for (u0,r0),(u1,r1) in zip(clean,clean[1:]):
            segment_lengths.append(quad(lambda t:math.hypot(r1-r0,(r0+t*(r1-r0))*(u1-u0)/R),
                                        0,1,epsabs=1e-11,epsrel=1e-12)[0])
        routes.append(dict(key=e['key'],face=e['layer'],source=source,target=target,
                           points_ur=clean,endpoints_ur=endpoints,
                           segment_arc_lengths=segment_lengths,length=sum(segment_lengths),
                           input_length=e['length']))
    devices=[dict(index=[r,c,0,j],radius=radius(r,c,j))
             for r in range(2) for c in range(12) for j in range(q)]
    return dict(schema='continuous-polar-routes-v2',code=code,p=p,q=q,M=12,
                minimum_qubit_pitch=1.,R=R,pitch=pitch,period=p*pitch,
                port_r_half=dr,port_u_half=du,devices=devices,orbits=routes,
                total_devices=p*len(devices),total_connections=p*len(routes),
                interpretation='A segment linearly interpolates u and r, then maps to x=r cos(u/R), y=r sin(u/R). Each orbit has all p rotations. Endpoints are (r,c,i,j). No intermediate face switches.',
                endpoint_exception='Only portions inside a fixed common endpoint region may merge for the corresponding pair of routes. Disjoint-endpoint routes are checked in full.',
                region='For center (u0,r0), |u-u0|<=port_u_half and |r-r0|<=port_r_half in reference unfolded coordinates, with periodic u.')

def stats(model,tier_faces):
    byface={}
    for e in model['orbits']:byface.setdefault(e['face'],[]).append(e)
    records=[];total=model['total_connections'];p=model['p']
    for t,faces in enumerate(tier_faces):
        edges=[e for f in faces if f is not None for e in byface[f]]
        n=len(edges)*p;s=sum(e['length'] for e in edges)*p
        upper=len(byface[faces[1]])*p if len(faces)==2 and faces[1] is not None else 0
        records.append(dict(tier=t,faces=faces,connections=n,mean_length=s/n,total_length=s,
                            TSV_total=2*t*n,HAL_Bonds=3*upper/n,
                            physical_face_switches=2*upper,
                            full_stack_bonds=(2*(t-1)*n+2*upper if t else 0)))
    T=len(records);L=sum(x['mean_length'] for x in records)/T
    tsv=sum(x['TSV_total'] for x in records);V=tsv/(total-records[0]['connections'])
    B=max(x['HAL_Bonds'] for x in records)
    counts=Counter((e['face'],tuple(ep)) for e in model['orbits'] for ep in (e['source'][:2]+[e['source'][3]],e['target'][:2]+[e['target'][3]]))
    return dict(tiers=T,faces=len(byface),Length=L,TSVs=V,Bonds=B,C_hw=1+(T-1)/16+(L-1)/36+V/12+B/16,
                total_edges=total,all_edge_mean=sum(e['length'] for e in model['orbits'])/len(model['orbits']),
                maximum_length=max(e['length'] for e in model['orbits']),TSV_total=tsv,
                physical_face_switches=sum(x['physical_face_switches'] for x in records),
                full_stack_bonds=sum(x['full_stack_bonds'] for x in records),
                max_face_endpoint_degree=max(counts.values()),tier_records=records)

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser();p.add_argument('geometry');p.add_argument('hardware');p.add_argument('output')
    a=p.parse_args();m=build(read(a.geometry),read(a.hardware));write(a.output,m)
    print('Explicit orbits',len(m['orbits']),'connections',m['total_connections'],
          'max length discrepancy',max(abs(e['length']-e['input_length']) for e in m['orbits']))
