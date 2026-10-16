#!/usr/bin/env python3
"""Build proMotif's promoter background for discovery mode.

For N random (or all) human protein-coding genes (Ensembl canonical transcript, main
chromosomes), take the promoter window -2500/+500 bp around the TSS (same as
proMotif's default), fetch the sequence (UCSC API, hg38), score every JASPAR CORE vertebrate matrix of a release on
both strands, keep each promoter's best-site log-odds score, and store the
empirical distribution (upper-tail quantiles + exponential tail) per matrix.

Matrices are converted exactly as proMotif's .pfm_to_logodds(): pseudocount 0.8
spread over a flat background, log2 odds vs 0.25.

Inputs : canon.tsv (Ensembl BioMart: gene id, name, chr, strand, tx start, tx end,
         transcript id; protein_coding, transcript_is_canonical), and the JASPAR
         TRANSFAC export of the release.
Output : promoter_background_JASPAR<release>_hg38.tsv
Usage  : python3 make_promoter_background.py canon.tsv JASPAR2026_..._transfac.txt 2026 [N|all]
"""
import sys, json, os, random, time, math, urllib.request
import numpy as np

canon, transfac, release = sys.argv[1], sys.argv[2], sys.argv[3]
N_ARG = sys.argv[4] if len(sys.argv) > 4 else "2000"   # number of genes, or "all"
UP, DOWN, SEED = 2500, 500, 1
CHR = {str(i) for i in range(1, 23)} | {"X", "Y"}

# ---- 1. sample genes ------------------------------------------------------
rows = [l.rstrip("\n").split("\t") for l in open(canon)][1:]
rows = [r for r in rows if r[2] in CHR]
seen, genes = set(), []
for r in rows:
    if r[0] not in seen:
        seen.add(r[0]); genes.append(r)
N = len(genes) if N_ARG == "all" else min(int(N_ARG), len(genes))
random.seed(SEED)
genes = random.sample(sorted(genes), N)

# ---- 2. fetch promoter sequences (UCSC REST API, hg38 = GRCh38, parallel) --
from concurrent.futures import ThreadPoolExecutor
cache = "bg_sequences.json"
seqs = json.load(open(cache)) if os.path.exists(cache) else {}
COMP = str.maketrans("ACGTN", "TGCAN")
def fetch(r):
    strand, s, e = int(r[3]), int(r[4]), int(r[5])
    tss = s if strand == 1 else e
    lo, hi = (tss - UP, tss + DOWN) if strand == 1 else (tss - DOWN, tss + UP)
    url = (f"https://api.genome.ucsc.edu/getData/sequence?genome=hg38;"
           f"chrom=chr{r[2]};start={max(0, lo - 1)};end={hi}")
    for attempt in range(5):
        try:
            dna = json.load(urllib.request.urlopen(url, timeout=60))["dna"].upper()
            return r[6], (dna if strand == 1 else dna.translate(COMP)[::-1])
        except Exception:
            time.sleep(2 * (attempt + 1))
    return r[6], ""
todo = [g for g in genes if g[6] not in seqs]
with ThreadPoolExecutor(8) as ex:
    for k, (tid, dna) in enumerate(ex.map(fetch, todo)):
        seqs[tid] = dna
        if k % 200 == 0: print(f"fetched {k}/{len(todo)}", flush=True)
        if k % 1000 == 999: json.dump(seqs, open(cache, "w"))     # resume-safe
json.dump(seqs, open(cache, "w"))
seq_list = [seqs[g[6]] for g in genes if len(seqs.get(g[6], "")) >= 100]
print("promoters:", len(seq_list))

# ---- 3. JASPAR matrices (same conversion as .pfm_to_logodds) ---------------
mots, cur = [], None
for l in open(transfac):
    l = l.rstrip("\n")
    if l.startswith("AC "): cur = {"id": l[3:].strip(), "rows": []}
    elif l.startswith("ID ") and cur: cur["name"] = l[3:].strip()
    elif cur and l[:2].isdigit(): cur["rows"].append([float(x) for x in l.split()[1:5]])
    elif l == "//" and cur: mots.append(cur); cur = None
def logodds(counts, pseudo=0.8):
    c = np.array(counts).T                      # 4 x L (A, C, G, T)
    ppm = (c + pseudo * 0.25) / (c.sum(0) + pseudo)
    return np.log2(ppm / 0.25)

# ---- 4. best-site score per promoter per matrix -----------------------------
code = {"A": 0, "C": 1, "G": 2, "T": 3}
enc = [np.array([code.get(ch, 4) for ch in s], dtype=np.int8) for s in seq_list]
mats = []
for m in mots:
    pwm = logodds(m["rows"])
    mats.append([np.vstack([P, np.zeros((1, pwm.shape[1]))]).astype(np.float32)
                 for P in (pwm, pwm[::-1, ::-1])])          # forward, reverse complement
# resumable: each block of promoters is saved to chunks/; finished blocks are skipped.
# Run several workers in parallel with WORKER=i NWORK=k; the last one assembles.
CHUNK = 500
os.makedirs("chunks", exist_ok=True)
WORKER, NWORK = int(os.environ.get("WORKER", 0)), int(os.environ.get("NWORK", 1))
t0 = time.time()
blocks = list(range(0, len(enc), CHUNK))
for bi, c0 in enumerate(blocks):
    out = f"chunks/best_{c0:06d}.npy"
    if bi % NWORK != WORKER or os.path.exists(out): continue
    part = enc[c0:c0 + CHUNK]
    blk = np.full((len(mots), len(part)), -np.inf, dtype=np.float32)
    cat = np.concatenate(part); starts = np.cumsum([0] + [len(e) for e in part])
    for mi, (fw, rc) in enumerate(mats):
        L = fw.shape[1]; n = len(cat) - L + 1
        for aug in (fw, rc):
            sc = np.zeros(n, dtype=np.float32)
            for j in range(L):
                sc += aug[cat[j:j + n], j]
            for gi in range(len(part)):
                a_, b_ = starts[gi], starts[gi + 1] - L + 1
                if b_ > a_:
                    v = sc[a_:b_].max()
                    if v > blk[mi, gi]: blk[mi, gi] = v
    np.save(out + ".tmp.npy", blk); os.replace(out + ".tmp.npy", out)
    print(f"block {c0} done  {time.time() - t0:.0f}s", flush=True)
missing = [c0 for c0 in blocks if not os.path.exists(f"chunks/best_{c0:06d}.npy")]
if missing:
    print(f"{len(missing)} blocks still missing"); sys.exit(0)
best = np.concatenate([np.load(f"chunks/best_{c0:06d}.npy") for c0 in blocks], axis=1)
best = best.astype(np.float64)

np.save(f"bg_best_JASPAR{release}.npy", best)          # raw maxima, for re-analysis

# ---- 5. empirical background per matrix -------------------------------------
# P(best site of a random promoter >= s) is read from upper-tail quantiles q_p
# (dense in the tail); beyond the strongest background promoter an exponential
# tail (peaks over threshold: mean excess over the top-5% quantile) is used.
# A Gumbel fit is NOT used for p-values: short motifs reach their maximum score
# in many promoters and are badly described by a smooth extreme-value curve.
# beta (Gumbel scale, moments) is kept only to shift the background for other
# window lengths (location + beta * log(L / L0)).
PROBS = [1.0, 0.5, 0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.075, 0.05, 0.04, 0.03,
         0.025, 0.02, 0.015, 0.01, 0.0075, 0.005, 0.004, 0.003, 0.0025, 0.002,
         0.0015, 0.001, 0.0005, 0.0002, 0.0001]
nprom = best.shape[1]
PROBS = [p for p in PROBS if p > 1.5 / nprom] + [1.0 / nprom]   # last = the floor, 1/N
with open(f"promoter_background_JASPAR{release}_hg38.tsv", "w") as f:
    f.write("\t".join(["matrix_id", "name", "n_bg", "mean", "sd", "beta", "tail_sigma"] +
                      [f"q_{p:g}" for p in PROBS]) + "\n")
    for mi, m in enumerate(mots):
        x = np.sort(best[mi][np.isfinite(best[mi])]); n = len(x)
        q = [x[min(n - 1, max(0, n - int(round(p * n))))] for p in PROBS]   # P(X >= q_p) = p
        u = q[PROBS.index(0.05)]; exc = x[x > u] - u
        sigma = max(exc.mean() if len(exc) else 0.0, 0.05 * x.std(ddof=1), 1e-3)
        beta = x.std(ddof=1) * math.sqrt(6) / math.pi
        f.write("\t".join([m["id"], m["name"], str(n), f"{x.mean():.4f}",
                           f"{x.std(ddof=1):.4f}", f"{beta:.4f}", f"{sigma:.4f}"] +
                          [f"{v:.4f}" for v in q]) + "\n")
print("done")
