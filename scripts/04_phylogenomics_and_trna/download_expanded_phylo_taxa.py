import os
import urllib.request
import json

ref_dir = "/mnt/d/Mitochondria/reference/complete_mitogenomes"
os.makedirs(ref_dir, exist_ok=True)

def download_gb(acc, out_name):
    out_path = os.path.join(ref_dir, out_name)
    if os.path.exists(out_path) and os.path.getsize(out_path) > 50000:
        print(f"[SKIP] {out_name} already exists ({os.path.getsize(out_path):,} bytes)")
        return out_path
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id={acc}&rettype=gbwithparts&retmode=text"
    print(f"[FETCH] Downloading {acc} -> {out_name} ...")
    req = urllib.request.Request(url, headers={"User-Agent": "AntigravityResearch/1.0"})
    data = urllib.request.urlopen(req, timeout=60).read()
    with open(out_path, "wb") as f:
        f.write(data)
    print(f"[DONE] Saved {out_name} ({len(data):,} bytes)")
    return out_path

# 1. Carya illinoinensis (Chen et al., 2024 Scientific Reports)
download_gb("PX519271.1", "Outgroup_Carya_illinoinensis_PX519271.1.gb")

# 2. Search for Fagus sylvatica complete mitochondrion
try:
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=Fagus+sylvatica[Organism]+AND+mitochondrion[Title]+AND+complete+genome[Title]&retmode=json"
    res = json.loads(urllib.request.urlopen(url, timeout=15).read())
    ids = res["esearchresult"]["idlist"]
    if ids:
        sum_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&id={ids[0]}&retmode=json"
        acc = json.loads(urllib.request.urlopen(sum_url, timeout=15).read())["result"][ids[0]]["accessionversion"]
        print(f"Found Fagus sylvatica mitogenome: {acc}")
        download_gb(acc, f"Outgroup_Fagus_sylvatica_{acc}.gb")
except Exception as e:
    print("Fagus sylvatica search/download failed:", e)

# 3. Arabidopsis thaliana complete mitochondrion (NC_037304.1)
download_gb("NC_037304.1", "Outgroup_Arabidopsis_thaliana_NC_037304.1.gb")
