"""Exactly optimize whole-face stacking with fixed continuous routes.

No route, code, pitch, clearance or edge set is changed. For each base face,
enumerate every pairing of the remaining faces (one dummy if needed). Put
the denser face on the lower surface and sort pairs by descending edge count.
These two operations independently minimize HAL Bonds and TSVs for a pairing.
"""
from pathlib import Path
import argparse
import hashlib
import json
import math
import time

HERE = Path(__file__).resolve().parent


def partitions(items):
    if not items:
        yield ()
    else:
        a = items[0]
        for i in range(1, len(items)):
            for rest in partitions(items[1:i] + items[i+1:]):
                yield ((a, items[i]),) + rest


def count_matchings(n):
    return math.prod(range(1, n, 2))


def cost(T, L, V, B):
    return 1 + ((T-1)/4 + max(0, L-1)/9 + V/3 + B/4)/4


def optimize(hw, source='unspecified'):
    counts = {int(k): v for k, v in hw['layer_counts'].items()}
    sums = {int(k): v for k, v in hw['layer_length_sums'].items()}
    ids = tuple(sorted(counts))
    total_n = sum(counts.values())
    total_length = sum(sums.values())
    T = 1 + math.ceil((len(ids)-1)/2)
    counts[-1] = 0
    sums[-1] = 0.

    def summarize(base, pairs, sort=True):
        # Ordering minimizes TSVs by the rearrangement inequality. Putting
        # the smaller face on top minimizes each pair's switch contribution.
        pairs = tuple(tuple(sorted(pair, key=lambda x: (-counts[x], x))) for pair in pairs)
        if sort:
            pairs = tuple(sorted(pairs, key=lambda pair: (-sum(counts[x] for x in pair), pair)))
        groups = ((base,),) + pairs
        records = []
        for tier, faces in enumerate(groups):
            n = sum(counts[x] for x in faces)
            length = sum(sums[x] for x in faces)
            upper_n = counts[faces[1]] if len(faces) == 2 else 0
            records.append(dict(tier=tier, faces=[None if x == -1 else x for x in faces],
                                edges=n, length_sum=length, mean_length=length/n,
                                TSV_total=2*tier*n, source_segments=3*upper_n,
                                physical_face_switches=2*upper_n,
                                full_stack_bonds=2*(tier-1)*n+2*upper_n if tier else 0))
        L = sum(r['mean_length'] for r in records)/T
        V = sum(r['TSV_total'] for r in records)/(total_n-counts[base])
        B = max(r['source_segments']/r['edges'] for r in records)
        out = dict(base_face=base, tier_faces=[r['faces'] for r in records], tiers=T,
                   Length=L, TSVs=V, Bonds=B, C_hw=cost(T,L,V,B),
                   TSV_total=sum(r['TSV_total'] for r in records),
                   physical_face_switches=sum(r['physical_face_switches'] for r in records),
                   full_stack_bonds=sum(r['full_stack_bonds'] for r in records),
                   all_edge_mean=total_length/total_n, all_edge_length_sum=total_length,
                   maximum_length=hw['maximum_length'], total_edges=total_n,
                   score_components=dict(tiers=(T-1)/16,length=max(0,L-1)/36,
                                         TSVs=V/12,Bonds=B/16), tier_records=records)
        assert sum(r['edges'] for r in records) == total_n
        assert math.isclose(sum(r['length_sum'] for r in records),total_length,rel_tol=1e-13)
        return out

    base_order = sorted(ids, key=lambda x: (-counts[x], x))
    remaining = base_order[1:]
    if len(remaining)%2:
        remaining.append(-1)
    baseline = summarize(base_order[0], tuple(tuple(remaining[i:i+2]) for i in range(0,len(remaining),2)),sort=False)

    # Cache pair summaries so enumeration uses only arithmetic, not dictionaries.
    all_ids=ids+(-1,)
    pc={}
    for i in all_ids:
        for j in all_ids:
            if i<j:
                n=counts[i]+counts[j]
                pc[(i,j)]=(n,(sums[i]+sums[j])/n,3*min(counts[i],counts[j])/n)
    best_cost = best_length = best_tsv = best_constrained = None
    feasible_constrained=0
    shortlist=[]
    weight_choices=(0.,.25,.5,1.,2.,4.,8.,16.)
    weighted_best={w:None for w in weight_choices}
    enumerated=0
    begin=time.perf_counter()
    for base in ids:
        remaining=tuple(x for x in ids if x!=base)
        if len(remaining)%2:
            remaining=(-1,)+remaining
        for pairs in partitions(remaining):
            d=[pc[tuple(sorted(pair))] for pair in pairs]
            L=(sums[base]/counts[base]+sum(x[1] for x in d))/T
            TSV_total=sum(2*t*n for t,n in enumerate(sorted((x[0] for x in d),reverse=True),1))
            V=TSV_total/(total_n-counts[base])
            B=max(x[2] for x in d)
            C=cost(T,L,V,B)
            item=(C,L,V,B,base,pairs)
            if best_cost is None or item[:4]<best_cost[:4]:
                best_cost=item
            if best_length is None or (L,C,V,B)<(best_length[1],best_length[0],best_length[2],best_length[3]):
                best_length=item
            if best_tsv is None or (V,C,L,B)<(best_tsv[2],best_tsv[0],best_tsv[1],best_tsv[3]):
                best_tsv=item
            full_stack_bonds=TSV_total-2*(total_n-counts[base])+sum(2*min(counts[i],counts[j]) for i,j in pairs)
            if TSV_total<=baseline['TSV_total'] and full_stack_bonds<=baseline['full_stack_bonds']:
                feasible_constrained+=1
                if best_constrained is None or item[:4]<best_constrained[:4]:
                    best_constrained=item
            for w in weight_choices:
                objective=w*(L-1)/36+V/12+B/16
                if weighted_best[w] is None or (objective,C)<weighted_best[w][:2]:
                    weighted_best[w]=(objective,C,item)
            enumerated+=1
    expected=len(ids)*count_matchings(len(ids)-1 if len(ids)%2 else len(ids))
    assert enumerated==expected,(enumerated,expected)
    def expand(item):
        result=summarize(item[4],item[5])
        assert math.isclose(result['C_hw'],item[0],rel_tol=1e-13)
        return result
    return dict(source=source,status='EXHAUSTIVE_FIXED_FACE_STACKING',
                continuous_arc_length_unit='minimum physical qubit pitch a0=1',
                number_of_nonempty_faces=len(ids),number_of_partitions=enumerated,
                elapsed_seconds=time.perf_counter()-begin,
                baseline=baseline,minimum_cost=expand(best_cost),minimum_Length=expand(best_length),
                minimum_TSVs=expand(best_tsv),
                no_more_physical_resources=dict(
                    constraints=dict(TSV_total_max=baseline['TSV_total'],full_stack_bonds_max=baseline['full_stack_bonds']),
                    feasible_partitions=feasible_constrained,
                    minimum_cost=expand(best_constrained) if best_constrained else None),
                sensitivity=[dict(length_term_weight=w,result=expand(weighted_best[w][2])) for w in weight_choices],
                scope='Whole existing faces moved intact. In-plane routes and all physical lengths are unchanged. No newly routed connections. Exact minimum only within this stacking family at fixed tier count.',
                physical_limitations='Inherited terminal, finite-width, vertical-via and device-region assumptions remain. No control/readout/grounding routing added.')
