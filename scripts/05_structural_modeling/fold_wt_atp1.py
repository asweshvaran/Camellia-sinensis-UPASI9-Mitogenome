import os
import torch
from Bio import SeqIO
from transformers import AutoTokenizer, EsmForProteinFolding

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
print(f"CUDA Available: {torch.cuda.is_available()} | Device: {torch.cuda.get_device_name(0) if torch.cuda.is_available() else 'CPU'}")

ref_faa = "/mnt/d/Mitochondria/reference/ref_mitochondrial_proteins.faa"
wt_seq = None
for r in SeqIO.parse(ref_faa, "fasta"):
    if r.id.startswith("atp1|"):
        wt_seq = str(r.seq).rstrip("*")
        break

assert wt_seq, "Reference atp1 not found!"
print(f"Loaded WT Reference ATP1: {len(wt_seq)} aa")

out_dir = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/gpu_structures"
os.makedirs(out_dir, exist_ok=True)
pdb_path = os.path.join(out_dir, "NC_043914_WT_atp1_predicted.pdb")

tokenizer = AutoTokenizer.from_pretrained("facebook/esmfold_v1")
model = EsmForProteinFolding.from_pretrained("facebook/esmfold_v1", torch_dtype=torch.bfloat16, low_cpu_mem_usage=True)
model = model.to(device)
model.eval()
model.trunk.set_chunk_size(64)

with torch.no_grad(), torch.amp.autocast("cuda", dtype=torch.bfloat16):
    inputs = tokenizer([wt_seq], return_tensors="pt", add_special_tokens=False)
    inputs = {k: v.to(device) for k, v in inputs.items()}
    output = model(**inputs)

    output_fp32 = {k: (v.float() if torch.is_tensor(v) and v.is_floating_point() else v) for k, v in output.items()}
    pdb_str = model.output_to_pdb(output_fp32)[0]
    with open(pdb_path, "w") as f:
        f.write(pdb_str)
    print(f"Saved WT ATP1 PDB cleanly to {pdb_path}")

    # Also fold recovered COX1 (527 aa) while model is warm!
    faa_path = "/mnt/d/Mitochondria/analysis/SRR9057080/annotation/predicted_proteins.faa"
    cox1_seq = None
    for r in SeqIO.parse(faa_path, "fasta"):
        if r.id == "NODE_15_length_1916_cov_19.127554_1":
            cox1_seq = str(r.seq).rstrip("*")
            break
    if cox1_seq:
        cox1_pdb = os.path.join(out_dir, "NODE_15_COX1_predicted.pdb")
        torch.cuda.empty_cache()
        with torch.no_grad(), torch.amp.autocast("cuda", dtype=torch.bfloat16):
            inputs = tokenizer([cox1_seq], return_tensors="pt", add_special_tokens=False)
            inputs = {k: v.to(device) for k, v in inputs.items()}
            output2 = model(**inputs)
        output2_fp32 = {k: (v.float() if torch.is_tensor(v) and v.is_floating_point() else v) for k, v in output2.items()}
        with open(cox1_pdb, "w") as f:
            f.write(model.output_to_pdb(output2_fp32)[0])
        print(f"Saved recovered COX1 PDB cleanly to {cox1_pdb}")
