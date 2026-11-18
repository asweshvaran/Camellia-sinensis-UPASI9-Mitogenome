import os
import csv
import urllib.request
from collections import Counter, defaultdict, deque
from Bio import SeqIO

# 1. Parse SPAdes GFA (assembly_graph_with_scaffolds.gfa) and GetOrganelle selected GFA
def parse_gfa(gfa_path):
    segs = {}
    links = []
    adj = defaultdict(set)
    paths = []
    with open(gfa_path) as f:
        for line in f:
            if line.startswith("S\t"):
                parts = line.strip().split("\t")
                sid = parts[1]
                seq = parts[2]
                slen = len(seq) if seq != "*" else 0
                cov = None
                for tag in parts[3:]:
                    if tag.startswith("LN:i:") and slen == 0:
                        slen = int(tag[5:])
                    if tag.startswith("dp:f:") or tag.startswith("KC:i:") or tag.startswith("RC:i:"):
                        cov = tag
                segs[sid] = slen
            elif line.startswith("L\t"):
                parts = line.strip().split("\t")
                u, v = parts[1], parts[3]
                links.append((u, v))
                adj[u].add(v)
                adj[v].add(u)
            elif line.startswith("P\t"):
                paths.append(line.strip().split("\t")[1])
    # Connected components
    visited = set()
    comps = []
    for s in segs:
        if s not in visited:
            q = deque([s])
            visited.add(s)
            comp = []
            while q:
                curr = q.popleft()
                comp.append(curr)
                for nb in adj[curr]:
                    if nb not in visited:
                        visited.add(nb)
                        q.append(nb)
            comps.append(comp)
    dead_ends = sum(1 for s in segs if len(adj[s]) <= 1)
    branching = sum(1 for s in segs if len(adj[s]) >= 3)
    isolated = sum(1 for c in comps if len(c) == 1)
    return {
        "nodes": len(segs),
        "total_bp": sum(segs.values()),
        "links": len(links),
        "components": len(comps),
        "largest_comp_nodes": max((len(c) for c in comps), default=0),
        "isolated_nodes": isolated,
        "dead_ends": dead_ends,
        "branching_hubs": branching,
        "paths": len(paths),
    }

gfa1 = "/mnt/d/Mitochondria/analysis/SRR9057080/assembly_mt/extended_spades/assembly_graph_with_scaffolds.gfa"
gfa2 = "/mnt/d/Mitochondria/analysis/SRR9057080/assembly_mt/embplant_mt.K105.contigs.graph1.selected_graph.gfa"
print("=== 1. ASSEMBLY GRAPH TOPOLOGY AUDIT ===")
print("  SPAdes full graph (assembly_graph_with_scaffolds.gfa):", parse_gfa(gfa1))
print("  GetOrganelle selected graph (K105 selected_graph.gfa):", parse_gfa(gfa2))

# 2. Fetch NCBI SRA RunInfo for all 28 accessions
print("\n=== 2. COHORT SRA METADATA VERIFICATION (28 ACCESSIONS) ===")
srr_file = "/mnt/d/Mitochondria/srr_list.txt"
with open(srr_file) as f:
    srrs = [l.strip() for l in f if l.strip()]

meta_out = "/mnt/d/Mitochondria/Q1_Investigation/Comparative_Genomics/cohort_28_sra_metadata.csv"
try:
    query = " OR ".join(srrs)
    url = f"https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/sra-db-be.cgi?rettype=runinfo&term={urllib.request.quote(query)}"
    req = urllib.request.urlopen(url, timeout=20)
    csv_text = req.read().decode("utf-8", errors="ignore")
    lines = [l for l in csv_text.splitlines() if l.strip()]
    if len(lines) > 1:
        with open(meta_out, "w", newline="") as out:
            out.write("\n".join(lines) + "\n")
        reader = csv.DictReader(lines)
        strat_cnt = Counter()
        src_cnt = Counter()
        bp_cnt = Counter()
        for r in reader:
            strat_cnt[r.get("LibraryStrategy", "?")] += 1
            src_cnt[r.get("LibrarySource", "?")] += 1
            bp_cnt[r.get("BioProject", "?")] += 1
        print(f"  Fetched metadata rows: {len(lines)-1}")
        print(f"  LibraryStrategy: {dict(strat_cnt)}")
        print(f"  LibrarySource  : {dict(src_cnt)}")
        print(f"  BioProjects    : {dict(bp_cnt)}")
    else:
        print("  Empty response from NCBI SRA runinfo.")
except Exception as e:
    print(f"  NCBI query failed: {e}")

# 3. Inspect Class C (NODE_2_ORF2) & Top Class D ORFs (NODE_10_ORF1, NODE_15_ORF1)
print("\n=== 3. SEQUENCE & BIOPHYSICAL CHECK OF KEY CANDIDATE ORFS ===")
faa_path = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/predicted_proteins.faa"
kd = {'I':4.5,'V':4.2,'L':3.8,'F':2.8,'C':2.5,'M':1.9,'A':1.8,'G':-0.4,'T':-0.7,'S':-0.8,'W':-0.9,'Y':-1.3,'P':-1.6,'H':-3.2,'E':-3.5,'Q':-3.5,'D':-3.5,'N':-3.5,'K':-3.9,'R':-4.5}
targets = [
    "NODE_2_length_3671_cov_47.563937_2",
    "NODE_10_length_2350_cov_28.737194_1",
    "NODE_15_length_1916_cov_19.127554_1",
    "NODE_1_length_7316_cov_129.055748_4",
]
for rec in SeqIO.parse(faa_path, "fasta"):
    if rec.id in targets:
        seq = str(rec.seq).rstrip("*")
        hydros = [kd.get(aa, 0.0) for aa in seq]
        # Simple 19-aa Kyte-Doolittle window TM prediction (>= 1.6 mean hydropathy)
        tm_windows = [sum(hydros[i:i+19])/19.0 for i in range(max(1, len(seq)-18))]
        max_tm = max(tm_windows) if tm_windows else 0
        mean_h = sum(hydros)/len(hydros) if hydros else 0
        print(f"  {rec.id} ({len(seq)} aa): mean_hydropathy={mean_h:.2f}, max_19aa_TM={max_tm:.2f} | start_desc={rec.description.split('#')[-1].strip()}")
        print(f"    N-term 60aa: {seq[:60]}")
