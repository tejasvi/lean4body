"""Re-verification of a certificate file written by bb.py, leaf by leaf, without the search:
    python verify_cert.py CERTFILE [nproc=10]
For every leaf (initial box + bisection path, rebuilt with the same float midpoints as bb.py):
  D leaf: the x-localization (bb.xenclose) keeps no interval, so the leaf contains no point of E;
  A leaf: every interval kept by the x-localization lies in one of the stored runs, and for every stored run the
          mean-value enclosure of tr S(y) on (leaf x run), with the STORED y, has a finite upper bound below the
          threshold (and not above the stored bound).
The structure (initial boxes = grid, complete trees) is checked by check_cert.py."""
import sys, os, math, time
import multiprocessing
FN = os.path.abspath(sys.argv[1])
NPROC = int(sys.argv[2]) if len(sys.argv) > 2 else 10
t = open(FN).readline().split()
MODE, NN, B0, PREC, THR = t[2], t[4:7], float(t[8]), t[10], t[12]
os.environ['PREC'] = PREC
os.environ['THRESH'] = THR
sys.argv = ['bb.py', MODE, NN[0], NN[1], NN[2], str(round(1 / B0)), '1e-5', '1']
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'work', 'rig'))
import bb, cert, jet
from flint import arb


def blocks():
    cur = None
    for line in open(FN):
        if line.startswith('#'):
            continue
        tok = line.split()
        if tok[0] == 'B':
            if cur is not None:
                yield cur
            v = [float(s) for s in tok[1:]]
            cur = ([(v[0], v[1]), (v[2], v[3]), (v[4], v[5])], [])
        elif tok[0] == 'D':
            cur[1].append(('' if tok[1] == '-' else tok[1], None))
        else:
            runs = [tuple(float(s) for s in r.split(',')) for r in tok[2:]]
            cur[1].append(('' if tok[1] == '-' else tok[1], runs))
    yield cur


def leafbox(bf0, path):
    bf = [list(lh) for lh in bf0]
    for i in range(0, len(path), 2):
        k, side = int(path[i]), path[i + 1]
        m = 0.5 * (bf[k][0] + bf[k][1])
        if side == 'L':
            bf[k][1] = m
        else:
            bf[k][0] = m
    return [tuple(lh) for lh in bf]


def tbound(bf, X, y):
    mids, rads, boxes = bb.box4(bf, X)
    yy = (arb(y[0]), arb(y[1]))
    tJ = cert.trS(bb.FUN(*jet.seeds(boxes)), yy)
    tm = cert.trS(bb.FUN(*[arb(m) for m in mids]), yy)
    v = jet.mv(tm, tJ, rads)
    return bb.up(v) if math.isfinite(float(v.mid())) and math.isfinite(float(v.rad())) else float('inf')


def check(block):
    bf0, leaves = block
    out = dict(nD=0, nA=0, nruns=0, bad=[], tmax=0.0)
    for path, runs in leaves:
        bf = leafbox(bf0, path)
        stats = dict(nP=0)
        kept, _ = bb.xenclose(bf, stats)
        if runs is None:
            out['nD'] += 1
            if kept is not None:
                out['bad'].append(('D leaf keeps an interval', bf0, path))
            continue
        out['nA'] += 1
        for lo, hi in (kept or []):
            if not any(r[0] <= lo and hi <= r[1] for r in runs):
                out['bad'].append(('kept interval outside the stored runs', bf0, path, (lo, hi)))
        for x0, x1, y0, y1, tu in runs:
            out['nruns'] += 1
            tb = tbound(bf, (x0, x1), (y0, y1))
            if not (tb < bb.THRESH and tb <= tu):
                out['bad'].append(('bound', bf0, path, tb, tu))
            out['tmax'] = max(out['tmax'], tb)
    return out


if __name__ == '__main__':
    t0 = time.time()
    tot = dict(nD=0, nA=0, nruns=0, bad=[], tmax=0.0, nB=0)
    with multiprocessing.get_context("fork").Pool(NPROC) as pool:
        for r in pool.imap_unordered(check, blocks(), chunksize=4):
            tot['nB'] += 1
            for k in ('nD', 'nA', 'nruns'):
                tot[k] += r[k]
            tot['bad'] += r['bad']
            tot['tmax'] = max(tot['tmax'], r['tmax'])
            if r['bad']:
                print("PROBLEM", r['bad'][:3], flush=True)
    print("%s %s: %d initial boxes, %d discarded leaves re-discarded, %d accepted leaves with %d runs re-verified; "
          "max bound %.12f < %s; problems %d; %.0f s" % (MODE, os.path.basename(FN), tot['nB'], tot['nD'], tot['nA'],
                                                        tot['nruns'], tot['tmax'], THR, len(tot['bad']), time.time() - t0))
