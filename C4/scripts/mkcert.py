"""Write the certificate modules `C4Cert/*.lean` from the hints found by `gen`.

    lake exe gen dir dir.txt && lake exe gen ch ch.txt
    python3 scripts/mkcert.py dir.txt ch.txt

Each line of a hint file is `n P C nodes clamps hint` for cell `n` (`P`, `C`: the numbers of
jets of `P` and of `tr S` evaluated by the replay, which predict its cost).  The cells are cut
into consecutive runs, one `decide +kernel` theorem per run (kept short, because the memory of a
kernel check grows with its length), and the runs into modules, which lake checks in parallel.
Each module chains its runs into a `Cover` of its range of cells, and `C4Cert/Cover.lean` chains
the modules into `dirCover : Cover dirCell 0 62720` and `chCover : Cover chCell 0 73728`.

It also writes `scripts/heavy_batches.txt` for `scripts/build_certs.sh`: the modules whose longest
theorem is predicted to take more than HEAVY_S seconds, in batches of at most BATCH_GB of kernel
memory, by the estimate 0.4 GB + 0.125 GB per second of the longest theorem.
"""
import os
import sys

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'C4Cert')
RUN_S = 6.0      # target predicted kernel time of one theorem, in seconds
MODULE_S = 100.0  # target predicted kernel time of one module
HEAVY_S = 8.0     # modules whose longest theorem is predicted to take longer are built first
BATCH_GB = 40.0   # memory budget of one batch of such modules


def cost(P, C):
    # measured: kernel seconds ~ 1.15 (0.07 P + 0.21 C) + 0.005 per cell
    return 1.15 * (0.07 * P + 0.21 * C) + 0.005


def load(path, N):
    rows = [l.split() for l in open(path)]
    assert len(rows) == N, (path, len(rows))
    for i, r in enumerate(rows):
        assert int(r[0]) == i, (path, i)
    return [(cost(int(r[1]), int(r[2])), r[5]) for r in rows]


def runs(cells):
    """consecutive runs `(a, b, cost)` of predicted cost about RUN_S"""
    out, a, c = [], 0, 0.0
    for n, (cn, _) in enumerate(cells):
        if n > a and c + cn > RUN_S:
            out.append((a, n, c))
            a, c = n, 0.0
        c += cn
    out.append((a, len(cells), c))
    return out


def modules(rs):
    out, cur, c = [], [], 0.0
    for r in rs:
        if cur and c + r[2] > MODULE_S:
            out.append(cur)
            cur, c = [], 0.0
        cur.append(r)
        c += r[2]
    out.append(cur)
    return out


def hint_lines(hs, indent='    '):
    lines, cur = [], indent
    for i, h in enumerate(hs):
        item = h + (',' if i + 1 < len(hs) else ']')
        if len(cur) > len(indent) and len(cur) + 1 + len(item) > 100:
            lines.append(cur)
            cur = indent
        cur += (' ' if len(cur) > len(indent) else '') + item
    lines.append(cur)
    return lines


def write_module(name, cell, cells, rs):
    a0, b0 = rs[0][0], rs[-1][1]
    L = ['import C4Check', '',
         f'/-! Cells `{a0} ≤ n < {b0}` of `{cell}`: `decide +kernel` checks that `checkBoxH` '
         'succeeds on each\nwith its hint. -/', '',
         'set_option Elab.async false', '',
         f'namespace C4.Cert.{name}', '']
    for k, (a, b, _) in enumerate(rs):
        L.append(f'theorem c{k} : allCells {cell} {a} {b} [')
        L.extend(hint_lines([cells[n][1] for n in range(a, b)]))
        L[-1] += ' = true := by'
        L.append('  decide +kernel')
        L.append('')
    chain = ' <|\n    '.join(f'(allCells_cover c{k}).trans' for k in range(len(rs) - 1))
    last = f'allCells_cover c{len(rs) - 1}'
    L.append(f'theorem cover : Cover {cell} {a0} {b0} :=')
    L.append('  ' + (chain + ' <|\n    ' + last if chain else last))
    L += ['', f'end C4.Cert.{name}', '']
    open(os.path.join(ROOT, name + '.lean'), 'w').write('\n'.join(L))
    return a0, b0


def peak_gb(max_s):
    return 0.4 + 0.125 * max_s


def heavy_batches(stats):
    """batches of the heavy modules, greedily by decreasing peak memory"""
    heavy = sorted((s for s in stats if s[2] > HEAVY_S), key=lambda s: -s[2])
    out, cur, mem = [], [], 0.0
    for s in heavy:
        if cur and mem + peak_gb(s[2]) > BATCH_GB:
            out.append(cur)
            cur, mem = [], 0.0
        cur.append(s[0])
        mem += peak_gb(s[2])
    if cur:
        out.append(cur)
    return out


def main():
    dir_path, ch_path = sys.argv[1], sys.argv[2]
    os.makedirs(ROOT, exist_ok=True)
    for f in os.listdir(ROOT):
        if f.endswith('.lean'):
            os.remove(os.path.join(ROOT, f))
    cover = {}
    stats = []
    for prefix, cell, path, N in [('Dir', 'dirCell', dir_path, 62720),
                                  ('Ch', 'chCell', ch_path, 73728)]:
        cells = load(path, N)
        ms = modules(runs(cells))
        names = []
        for i, rs in enumerate(ms):
            name = f'{prefix}{i:03d}'
            a0, b0 = write_module(name, cell, cells, rs)
            names.append((name, a0, b0))
            stats.append((name, sum(r[2] for r in rs), max(r[2] for r in rs), len(rs)))
        assert names[0][1] == 0 and names[-1][2] == N
        assert all(names[i][2] == names[i + 1][1] for i in range(len(names) - 1))
        cover[cell] = (names, N)
    L = [f'import C4Cert.{n}' for c in ('dirCell', 'chCell') for n, _, _ in cover[c][0]]
    L += ['', '/-!', '# The two covers', '',
          'Every cell of the two grids of `C4Check/Grid.lean` has a hint on which `checkBoxH`',
          'succeeds.  Each module `C4Cert.DirNNN`, `C4Cert.ChNNN` checks a range of cells by',
          '`decide +kernel` (written by `scripts/mkcert.py` from the hints found by `gen`).', '-/',
          '', 'namespace C4', '']
    for cell, cname in (('dirCell', 'dirCover'), ('chCell', 'chCover')):
        names, N = cover[cell]
        L.append(f'theorem {cname} : Cover {cell} 0 {N} :=')
        L.append('  ' + ' <|\n  '.join(f'Cert.{n}.cover.trans' for n, _, _ in names[:-1])
                 + ' <|\n  ' + f'Cert.{names[-1][0]}.cover')
        L.append('')
    L += ['end C4', '']
    open(os.path.join(ROOT, 'Cover.lean'), 'w').write('\n'.join(L))
    tot = sum(s[1] for s in stats)
    print(f'{len(stats)} modules, {sum(s[3] for s in stats)} theorems, predicted {tot:.0f} s')
    big = sorted(stats, key=lambda s: -s[2])
    print('largest theorems (module, module s, max theorem s):',
          [(s[0], round(s[1]), round(s[2], 1)) for s in big[:12]])
    peaks = sorted((peak_gb(s[2]) for s in stats), reverse=True)
    print(f'sum of the 32 largest module peaks: {sum(peaks[:32]):.1f} GB')
    batches = heavy_batches(stats)
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'heavy_batches.txt')
    open(path, 'w').write(''.join(' '.join(b) + '\n' for b in batches))
    print('heavy batches:', [len(b) for b in batches])


if __name__ == '__main__':
    main()
