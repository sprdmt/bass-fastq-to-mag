#!/usr/bin/env python3
"""Randomised comparison sweep behind the numbers in the GraphBin2 validation report.

Needs python-igraph and graphbin2 1.3.3 (for the previously used patch and unmodified upstream).
Usage: python tests/graphbin2_fast/sweep_graphbin2.py [--n 600] [--out sweep.json]
Seeds: tie-free 10000.., tie-heavy 10000.. (tie_free=False), upstream 20000..
"""
import argparse
import json
import statistics
import sys
import tempfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gb2_testlib as T  # noqa: E402


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=600)
    ap.add_argument("--n-upstream", type=int, default=400)
    ap.add_argument("--out", default="sweep_graphbin2.json")
    a = ap.parse_args()
    fast, ref = T.load_fast(), T.load_reference()
    ref_total, patch, up = T.load_reference(total_order=True), T.load_patch(), T.load_upstream()
    out = {}
    for mode in ("tie_free", "ties"):
        rows = []
        for seed in range(a.n):
            inst = T.random_instance(10000 + seed, tie_free=(mode == "tie_free"))
            with tempfile.TemporaryDirectory() as d:
                inst.write(d)
                rf = T.run_fast(fast, d, f"{d}/f")
                rr = T.run_module(ref, d, f"{d}/r")
                rt = T.run_module(ref_total, d, f"{d}/t")
                rp = T.run_module(patch, d, f"{d}/p") if patch else None
                rows.append(dict(eq_ref=rf == rr, eq_ref_total=rf == rt,
                                 agree_ref=T.agreement(rf, rr),
                                 eq_patch=(rp == rr) if rp is not None else None,
                                 agree_patch=T.agreement(rp, rr) if rp is not None else None,
                                 changed=rf != T.initial_result(inst)))
        out[mode] = dict(
            instances=a.n, fast_eq_reference=sum(r["eq_ref"] for r in rows),
            fast_eq_reference_total_order=sum(r["eq_ref_total"] for r in rows),
            mean_agreement_fast_reference=statistics.mean(r["agree_ref"] for r in rows),
            min_agreement_fast_reference=min(r["agree_ref"] for r in rows),
            patch_eq_reference=sum(1 for r in rows if r["eq_patch"]),
            mean_agreement_patch_reference=statistics.mean(
                r["agree_patch"] for r in rows if r["agree_patch"] is not None) if patch else None,
            instances_changed_by_graphbin2=sum(r["changed"] for r in rows))
    if up is not None:
        cnt = Counter()
        for seed in range(a.n_upstream):
            inst = T.random_instance(20000 + seed, tie_free=True)
            with tempfile.TemporaryDirectory() as d:
                inst.write(d)
                rr = T.run_module(ref, d, f"{d}/r")
                try:
                    ru = T.run_module(up, d, f"{d}/u")
                    cnt["ok_equal" if ru == rr else "ok_different"] += 1
                except BaseException as e:  # noqa: BLE001 - classify upstream failures
                    cnt[type(e).__name__] += 1
        out["upstream_unfixed_tie_free"] = dict(instances=a.n_upstream, outcomes=dict(cnt))
    Path(a.out).write_text(json.dumps(out, indent=1))
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
