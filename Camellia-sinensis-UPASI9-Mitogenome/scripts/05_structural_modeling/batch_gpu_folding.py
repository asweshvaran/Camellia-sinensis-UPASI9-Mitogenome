import os
import torch
from Bio import SeqIO
from transformers import AutoTokenizer, EsmForProteinFolding

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"CUDA Available: {torch.cuda.is_available()} | Device: {torch.cuda.get_device_name(0)}")

faa_path = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/predicted_proteins.faa"
out_dir = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/gpu_structures"
os.makedirs(out_dir, exist_ok=True)

# Top uncharacterized novel ORFs identified in our earlier discovery step
target_novel_ids = [
    "NODE_10_length_2350_cov_28.737194_1",  # 564 aa
    "NODE_15_length_1916_cov_19.127554_1",  # 528 aa
    "NODE_9_length_2846_cov_208.218533_1",  # 505 aa
]

records_map = {}
for r in SeqIO.parse(faa_path, "fasta"):
    for tid in target_novel_ids:
        if r.id.startswith(tid):
            records_map[tid] = str(r.seq).rstrip("*")

print(f"Found {len(records_map)} novel targets to fold on GPU.")

tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", low_cpu_mem_usage=True)
model = model.to(device)
model.eval()

for tid, seq in records_map.items():
    pdb_path = os.path.join(out_dir, f"{tid}_predicted.pdb")
    if os.path.exists(pdb_path):
        print(f"Skipping {tid}, already modeled.")
        continue
    
    print(f"Folding {tid} ({len(seq)} aa) on RTX 5070 Ti...")
    with torch.no_grad():
        inputs = tokenizer([seq], return_tensors="pt", add_special_tokens=False)
        inputs = {k: v.to(device) for k, v in inputs.items()}
        output = model(**inputs)
        
    pdb_str = model.output_to_pdb(output)[0]
    with open(pdb_path, "w") as f:
        f.write(pdb_str)
    print(f"Saved: {pdb_path}")

print("Batch GPU modeling finished successfully!")
