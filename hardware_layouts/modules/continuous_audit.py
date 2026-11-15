"""Independent complete periodic curve audit for the q3 final layouts.

The routes are rebuilt directly from geometry and physical edge records.
No sampled polygonal approximation to polar curves is used. Only COMMON
endpoint rectangles are removed for a pair of paths. Devices are checked
against complete routes, except the route's two actual endpoint devices.
"""
from pathlib import Path
import ast,math,json,itertools,time,hashlib,argparse

HERE=Path(__file__).resolve().parent

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def write(p,x):Path(p).write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n',encoding='utf8')

def point_distance(p,a,b):
    dx=b[0]-a[0];dy=b[1]-a[1];den=dx*dx+dy*dy
    t=max(0.,min(1.,((p[0]-a[0])*dx+(p[1]-a[1])*dy)/den)) if den else 0.
    q=(a[0]+t*dx,a[1]+t*dy)
    return math.hypot(p[0]-q[0],p[1]-q[1]),p,q

def segment_distance(a,b,c,d):
    u=(b[0]-a[0],b[1]-a[1]);v=(d[0]-c[0],d[1]-c[1]);z=(c[0]-a[0],c[1]-a[1])
    cross=u[0]*v[1]-u[1]*v[0]
    if abs(cross)>1e-16:
        t=(z[0]*v[1]-z[1]*v[0])/cross;s=(z[0]*u[1]-z[1]*u[0])/cross
        if -1e-13<=t<=1+1e-13 and -1e-13<=s<=1+1e-13:
            q=(a[0]+t*u[0],a[1]+t*u[1]);return 0.,q,q
    return min((point_distance(a,c,d),point_distance(b,c,d),point_distance(c,a,b),point_distance(d,a,b)),key=lambda z:z[0])

def clip(points,ports,dr,du):
    result=[]
    for a,b in zip(points,points[1:]):
        v=(b[0]-a[0],b[1]-a[1]);ts=[0.,1.]
        for x,r in ports:
            for dim,lo,hi in ((0,x-du,x+du),(1,r-dr,r+dr)):
                if abs(v[dim])>1e-15:
                    for z in (lo,hi):
                        t=(z-a[dim])/v[dim]
                        if 0<t<1:ts.append(t)
        ts=sorted(set(ts))
        for t0,t1 in zip(ts,ts[1:]):
            m=(a[0]+.5*(t0+t1)*v[0],a[1]+.5*(t0+t1)*v[1])
            if any(abs(m[0]-x)<du+1e-12 and abs(m[1]-r)<dr+1e-12 for x,r in ports):continue
            aa=(a[0]+t0*v[0],a[1]+t0*v[1]);bb=(a[0]+t1*v[0],a[1]+t1*v[1])
            if math.dist(aa,bb)>1e-13:result.append((aa,bb))
    return result

def rectangle_distance(a,b,x,r,dr,du):
    lo=(x-du,r-dr);hi=(x+du,r+dr)
    for q in (a,b):
        if lo[0]<=q[0]<=hi[0] and lo[1]<=q[1]<=hi[1]:return 0.,q,q
    corners=[lo,(hi[0],lo[1]),hi,(lo[0],hi[1])]
    return min((segment_distance(a,b,corners[i],corners[(i+1)%4]) for i in range(4)),key=lambda z:z[0])

def local_polar_factor(R,rmin,threshold):
    if threshold<=0:return min(1.,rmin/R)
    t=math.asin(min(1.,threshold/(2*rmin)))
    return min(1.,rmin/R*math.sin(t)/t)

def bound(d,rmin,R):
    if not math.isfinite(d):return math.inf
    v=min(rmin,d*min(1.,rmin/R))
    for _ in range(8):v=min(rmin,d*local_polar_factor(R,rmin,v))
    return max(0.,v-1e-10)

def reconstruct(g,h):
    cfg=g['parameters'];q=g['q'];p=g['p'];R=g.get('inner_radius',g.get('R'));pitch=cfg['pitch']
    assert abs(p*pitch-2*math.pi*R)<1e-8
    dr=cfg['terminal_r_half'];du=cfg['terminal_u_half'];lane=cfg['lane_spacing']
    rho=g.get('rho',[0,2,4,6,8,10,11,9,7,5,3,1]);ori=g.get('orientation',[1]*6+[-1]*6)
    low={};high={};r=R
    widths={ast.literal_eval(key):w for key,w in g['local_widths'].items()}
    if 'low_radii' in g:
        for key,w in widths.items():
            low[key]=g['low_radii'][str(key)];high[key]=g.get('high_radii',{}).get(str(key),low[key]+w)
    else:
        order=g.get('physical_module_channel_order',[list(range(q)) for _ in range(12)])
        for c in sorted(range(12),key=lambda c:rho[c]):
            for j in order[c]:low[c,j]=r;high[c,j]=r+widths[c,j];r=high[c,j]+1.
    def radius(rail,c,j):
        if rail==0:return high[c,j] if ori[c]>0 else low[c,j]
        c=(c+1)%12;return low[c,j] if ori[c]>0 else high[c,j]
    actual={}
    for e in h['edge_orbits']:
        if e['key'][0]=='v':
            _,c,j,k=e['key'];actual.setdefault((c,j),set()).add(k)
    curves=[];groups={ast.literal_eval(key):ds for key,ds in g['groups'].items()}
    for e in h['edge_orbits']:
        kind,c,j,k=e['key'];a=radius(*e['source']);b=radius(*e['target']);sign=1 if b>a else -1
        if kind=='h':
            ds=next(ds for ds in groups[c,j] if k in ds)
            eta=cfg.get('pair_offset',.25)
            alpha=0. if len(ds)==1 else (-eta if k==min(ds) else eta)
            alpha=g.get('horizontal_alpha_by_key',{}).get(str(tuple(e['key'])),alpha)
            pts=[(0,a),(alpha*pitch,a),(alpha*pitch,a+sign*dr),
                 ((k-alpha)*pitch,b-sign*dr),((k-alpha)*pitch,b),(k*pitch,b)]
            ports=[(0,a),(k*pitch,b)]
        else:
            targets=actual[c,j] if g.get('vertical_rank_mode','all_targets')=='actual_targets' else range(q)
            order=sorted(targets,key=lambda k:abs(radius(1,c,k)-a));rank=order.index(k)
            center=g.get('vertical_offset_center',1.5)
            theta_sign=g.get('vertical_sign_override',{}).get(str(tuple(e['key'])),g.get('radial_angular_sign_by_source_channel',[(-1)**j for j in range(q)])[j])
            explicit=g.get('vertical_route_parameters',{}).get(str(tuple(e['key'])),{})
            offset=explicit.get('radial_offset',sign*(center-rank)*lane)
            u=explicit.get('angular_offset',theta_sign*(du+(rank+1)*lane));rp=a+offset
            pts=[(0,a),(0,rp),(u,rp),(u,b),(0,b)];ports=[(0,a),(0,b)]
        clean=[pts[0]]
        for z in pts[1:]:
            if math.dist(z,clean[-1])>1e-13:clean.append(z)
        curves.append(dict(key=e['key'],face=e['layer'],points=clean,ports=ports,
             umin=min(x[0] for x in clean),umax=max(x[0] for x in clean),
             rmin=min(x[1] for x in clean),rmax=max(x[1] for x in clean)))
    radii=sorted(list(low.values())+list(high.values()))
    return dict(R=R,pitch=pitch,q=q,p=p,dr=dr,du=du,low=low,high=high,radii=radii,curves=curves)

def audit(gpath,hpath):
    begin=time.perf_counter();g=read(gpath);h=read(hpath);data=reconstruct(g,h)
    R=data['R'];pitch=data['pitch'];dr=data['dr'];du=data['du'];p=data['p'];cs=data['curves'];radii=data['radii']
    stat={family:dict(orbit_pairs=0,periodic_pair_classes=0,segment_pair_checks=0,
          radial_exclusions=0,radial_lower=math.inf,unfolded_lower=math.inf,witness=None) for family in ('HH','VV','HV')}
    for i,a in enumerate(cs):
        for j in range(i,len(cs)):
            b=cs[j]
            if a['face']!=b['face']:continue
            family='HH' if a['key'][0]==b['key'][0]=='h' else 'VV' if a['key'][0]==b['key'][0]=='v' else 'HV'
            s=stat[family];s['orbit_pairs']+=1
            rgap=max(a['rmin'],b['rmin'])-min(a['rmax'],b['rmax'])
            if rgap>.100001:
                s['radial_exclusions']+=1;s['radial_lower']=min(s['radial_lower'],rgap);continue
            nlo=math.floor((a['umin']-b['umax'])/pitch)-1;nhi=math.ceil((a['umax']-b['umin'])/pitch)+1
            assert max(abs(nlo),abs(nhi))<p/2-2
            for n in range(nlo,nhi+1):
                if i==j and n<=0:continue
                bp=[(u+n*pitch,r) for u,r in b['points']];be=[(u+n*pitch,r) for u,r in b['ports']]
                common=[x for x in a['ports'] if any(math.dist(x,y)<1e-10 for y in be)]
                ac=clip(a['points'],common,dr,du);bc=clip(bp,common,dr,du)
                s['periodic_pair_classes']+=1
                for ab,cd in itertools.product(ac,bc):
                    d,x,y=segment_distance(*ab,*cd);s['segment_pair_checks']+=1
                    if d<s['unfolded_lower']:
                        s['unfolded_lower']=d;s['witness']=dict(edge1=a['key'],edge2=b['key'],station2=n,face=a['face'],shared_ports=common,segment1=ab,segment2=cd,point1=x,point2=y)
    # The angular windows omitted by the integer interval have at least one
    # reference pitch of angular separation, including their periodic copies.
    rmin=min(c['rmin'] for c in cs);far_wire=2*rmin*math.sin(pitch/(2*R))
    for family,s in stat.items():
        s['continuous_lower']=min(bound(s['unfolded_lower'],rmin,R),s['radial_lower'],far_wire)
        s['status']='PASS' if s['continuous_lower']>.1 else 'FAIL'
        for key in ('unfolded_lower','radial_lower'):
            if not math.isfinite(s[key]):s[key]=None
    devs={family:dict(segment_checks=0,radial_exclusions=0,radial_lower=math.inf,unfolded_lower=math.inf,witness=None) for family in ('H','V')}
    for a in cs:
        s=devs[a['key'][0].upper()]
        nlo=math.floor((a['umin']-du)/pitch)-1;nhi=math.ceil((a['umax']+du)/pitch)+1
        assert max(abs(nlo),abs(nhi))<p/2-2
        for r in radii:
            rgap=max(a['rmin'],r-dr)-min(a['rmax'],r+dr)
            if rgap>.100001:
                s['radial_exclusions']+=1;s['radial_lower']=min(s['radial_lower'],rgap);continue
            for n in range(nlo,nhi+1):
                x=n*pitch
                if any(math.dist((x,r),z)<1e-10 for z in a['ports']):continue
                for ab in zip(a['points'],a['points'][1:]):
                    d,xx,yy=rectangle_distance(*ab,x,r,dr,du);s['segment_checks']+=1
                    if d<s['unfolded_lower']:
                        s['unfolded_lower']=d;s['witness']=dict(edge=a['key'],device=(x,r),segment=ab,point1=xx,point2=yy)
    far_device=2*(min(radii)-dr)*math.sin(pitch/(2*R))
    for family,s in devs.items():
        s['continuous_lower']=min(bound(s['unfolded_lower'],min(radii)-dr,R),s['radial_lower'],far_device)
        s['status']='PASS' if s['continuous_lower']>.05 else 'FAIL'
        for key in ('unfolded_lower','radial_lower'):
            if not math.isfinite(s[key]):s[key]=None
    minimum_pitch=min(min(b-a for a,b in zip(radii,radii[1:])),2*min(radii)*math.sin(math.pi/p))
    return dict(status='PASS' if all(z['status']=='PASS' for z in list(stat.values())+list(devs.values())) and abs(minimum_pitch-1)<1e-10 else 'FAIL',
       geometry_sha256=hashlib.sha256(Path(gpath).read_bytes()).hexdigest(),hardware_sha256=hashlib.sha256(Path(hpath).read_bytes()).hexdigest(),
       p=p,q=data['q'],orbit_count=len(cs),physical_edges=p*len(cs),device_centers=p*len(radii),
       minimum_qubit_center_distance=minimum_pitch,wire=stat,device=devs,
       minimum_wire_lower_bound=min(z['continuous_lower'] for z in stat.values()),
       minimum_device_region_lower_bound=min(z['continuous_lower'] for z in devs.values()),
       far_station_bounds=dict(wire=far_wire,device=far_device),
       boundary_policy='Only common endpoint regions clipped for pairwise wire checks. Complete routes tested against every nonendpoint device region.',
       method='Independent polar curve reconstruction. Exact unfolded segment distances and local hypothetical-violation polar contraction bound. All relevant periodic offsets plus explicit angular/radial bounds for excluded classes.',
       numerical_scope='Double precision exhaustive finite geometric cases and analytical bounds, not interval-arithmetic formal proof. Not a manufacturing, control/readout or electromagnetic validation.',seconds=time.perf_counter()-begin)

def compare_model(gpath,hpath,model):
    ds=reconstruct(read(gpath),read(hpath));assert len(ds['curves'])==len(model['orbits'])
    worst=0.
    for c,e in zip(ds['curves'],model['orbits']):
        assert c['key']==e['key'] and c['face']==e['face']
        assert len(c['points'])==len(e['points_ur']),(c['key'],c['points'],e['points_ur'])
        worst=max(worst,max(math.dist(a,b) for a,b in zip(c['points'],e['points_ur'])))
    return dict(status='PASS' if worst<1e-11 else 'FAIL',orbits=len(ds['curves']),maximum_point_discrepancy=worst)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('geometry',type=Path);p.add_argument('hardware',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
    x=audit(a.geometry,a.hardware);write(a.output,x)
    print(json.dumps({k:x[k] for k in ('status','minimum_wire_lower_bound','minimum_device_region_lower_bound','seconds')},indent=2))
