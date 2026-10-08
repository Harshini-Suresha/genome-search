"""Dynamic-programming alignment: Smith-Waterman (local) and Needleman-Wunsch (global), linear gap penalty."""
import numpy as np

def align(a, b, mode='local', match=2, mismatch=-1, gap=-2):
    n, m, loc = len(a), len(b), mode == 'local'
    H = [[0] * (m + 1) for _ in range(n + 1)]; T = [[0] * (m + 1) for _ in range(n + 1)]
    if not loc:
        for i in range(1, n + 1): H[i][0] = i * gap; T[i][0] = 2
        for j in range(1, m + 1): H[0][j] = j * gap; T[0][j] = 3
    best, bi, bj = 0, 0, 0
    for i in range(1, n + 1):
        ai, Hi, Hp, Ti = a[i-1], H[i], H[i-1], T[i]
        for j in range(1, m + 1):
            d = Hp[j-1] + (match if ai == b[j-1] else mismatch); u = Hp[j] + gap; l = Hi[j-1] + gap
            v = max(d, u, l); t = 1 if v == d else 2 if v == u else 3
            if loc and v <= 0: v, t = 0, 0
            Hi[j] = v; Ti[j] = t
            if loc and v > best: best, bi, bj = v, i, j
    if not loc: best, bi, bj = H[n][m], n, m
    i, j, ra, rb, rm, path = bi, bj, [], [], [], []
    while i > 0 or j > 0:
        if loc and H[i][j] == 0: break
        t = T[i][j]; path.append((i, j))
        if t == 1: ra.append(a[i-1]); rb.append(b[j-1]); rm.append('|' if a[i-1] == b[j-1] else 'x'); i -= 1; j -= 1
        elif t == 2: ra.append(a[i-1]); rb.append('-'); rm.append(' '); i -= 1
        elif t == 3: ra.append('-'); rb.append(b[j-1]); rm.append(' '); j -= 1
        else: break
    ra, rb, rm = ''.join(ra[::-1]), ''.join(rb[::-1]), ''.join(rm[::-1])
    L = max(1, len(rm))
    return dict(score=best, a=ra, b=rb, mid=rm, H=np.array(H), path=path[::-1], end=(bi, bj), start=(i, j),
                identity=rm.count('|') / L, gaps=ra.count('-') + rb.count('-'), length=len(rm))

def score_only(a, b, match=2, mismatch=-1, gap=-2):
    """Two-row local alignment score: O(mn) time, O(n) memory."""
    prev = [0] * (len(b) + 1); best = 0
    for ai in a:
        cur = [0]; c0 = 0
        for j, bj in enumerate(b, 1):
            v = max(0, prev[j-1] + (match if ai == bj else mismatch), prev[j] + gap, c0 + gap)
            cur.append(v); c0 = v
            if v > best: best = v
        prev = cur
    return best
