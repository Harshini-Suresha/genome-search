"""Sequence utilities, k-mer index, trie and reference scan used by the Genome Search workbench."""
import random, time, tracemalloc
from collections import Counter

def clean(s):
    return ''.join(c for c in s.upper() if c.isalpha())

def kind_of(s):
    return 'DNA' if set(s) <= set('ACGTN') else 'PROTEIN'

def gc(s):
    return (s.count('G') + s.count('C')) / len(s) if s else 0.0

def parse_fasta(text):
    out, name, buf = {}, None, []
    def flush():
        if name is not None and buf:
            key, i = name, 2
            while key in out: key, i = f"{name}_{i}", i + 1
            out[key] = clean(''.join(buf))
    for line in text.splitlines():
        line = line.strip()
        if not line: continue
        if line.startswith('>'):
            flush(); name = (line[1:].split() or [f'seq{len(out)+1}'])[0]; buf = []
        else:
            if name is None: name = 'seq1'
            buf.append(line)
    flush()
    return {k: v for k, v in out.items() if v}

def mutate(s, rate, rng, indel=0.01):
    out = []
    for c in s:
        r = rng.random()
        if r < indel / 2: continue
        if r < indel: out.append(c); out.append(rng.choice('ACGT')); continue
        out.append(rng.choice([x for x in 'ACGT' if x != c]) if r < indel + rate else c)
    return ''.join(out)

def synthetic(n=12, length=2000, gc_target=.5, motif='GATTACAGATTACA', related=.4, rate=.6, seed=1):
    rng = random.Random(seed); w = [(1-gc_target)/2, gc_target/2, gc_target/2, (1-gc_target)/2]
    base = ''.join(rng.choices('ACGT', w, k=length)); out = {}
    for i in range(n):
        s = base if i == 0 else (mutate(base, .08, rng) if rng.random() < related else ''.join(rng.choices('ACGT', w, k=length)))
        if motif and rng.random() < rate and len(s) > len(motif) + 1:
            p = rng.randrange(0, len(s) - len(motif)); s = s[:p] + motif + s[p+len(motif):]
        out[f'seq_{i+1:03d}'] = s
    return out

def naive_find(seq, q):
    out, i = [], seq.find(q)
    while i != -1: out.append(i); i = seq.find(q, i + 1)
    return out

def hamming(a, b):
    return sum(x != y for x, y in zip(a, b))

class KmerIndex:
    """Hash table from every k-mer to (sequence, position) pairs."""
    def __init__(self, k):
        self.k, self.table, self.names, self.seqs = k, {}, [], []
    def add(self, name, seq):
        sid = len(self.names); self.names.append(name); self.seqs.append(seq); k, t = self.k, self.table
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer in t: t[kmer].append((sid, i))
            else: t[kmer] = [(sid, i)]
    def exact(self, q):
        if len(q) < self.k: return []
        return [(self.names[s], p) for s, p in self.table.get(q[:self.k], ()) if self.seqs[s].startswith(q, p)]
    def approx(self, q, d):
        """Seed-and-vote: every query k-mer votes for a start position; candidates are verified by Hamming distance."""
        k, votes = self.k, Counter()
        for j in range(len(q) - k + 1):
            for s, p in self.table.get(q[j:j+k], ()):
                votes[(s, p - j)] += 1
        rows = []
        for (s, st), v in votes.items():
            seq = self.seqs[s]
            if st < 0 or st + len(q) > len(seq): continue
            mm = hamming(q, seq[st:st+len(q)])
            if mm <= d: rows.append((self.names[s], st, mm, v))
        return sorted(rows, key=lambda r: (r[2], r[0], r[1]))
    def guaranteed(self, qlen, d):
        """With d mismatches at most d*k k-mers are broken, so one survives if qlen-k+1 > d*k."""
        return qlen - self.k + 1 > d * self.k
    def stats(self):
        n = sum(len(v) for v in self.table.values())
        return dict(distinct=len(self.table), entries=n, avg_bucket=n / max(1, len(self.table)), max_bucket=max((len(v) for v in self.table.values()), default=0))
    def bucket_sizes(self):
        return Counter(len(v) for v in self.table.values())

def build_index(seqs, k):
    tracemalloc.start(); t = time.perf_counter(); idx = KmerIndex(k)
    for n, s in seqs.items(): idx.add(n, s)
    secs = time.perf_counter() - t; _, peak = tracemalloc.get_traced_memory(); tracemalloc.stop()
    return idx, secs, peak / 1e6

class Trie:
    def __init__(self): self.root, self.nodes, self.words = {}, 1, 0
    def insert(self, word, payload):
        node = self.root
        for ch in word:
            nxt = node.get(ch)
            if nxt is None: nxt = node[ch] = {}; self.nodes += 1
            node = nxt
        node.setdefault('$', []).append(payload); self.words += 1
    def prefix(self, p, limit=500):
        node = self.root
        for ch in p:
            node = node.get(ch)
            if node is None: return []
        out, stack = [], [node]
        while stack and len(out) < limit:
            n = stack.pop()
            if '$' in n: out.extend(n['$'])
            stack.extend(v for k, v in n.items() if k != '$')
        return out[:limit]
