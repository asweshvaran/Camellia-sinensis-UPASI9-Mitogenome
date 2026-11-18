#!/bin/bash
set -e
mkdir -p /mnt/d/Mitochondria/sra_data
cd /mnt/d/Mitochondria/sra_data

LOG="/mnt/d/Mitochondria/download_sra.log"
echo "[$(date)] SRA download pipeline started." >> "$LOG"

for srr in $(cat /mnt/d/Mitochondria/srr_list.txt); do
    echo "[$(date)] Downloading $srr..." >> "$LOG"
    if [ ! -f "$srr/$srr.sra" ]; then
        prefetch --max-size 100G --output-directory . "$srr" >> "$LOG" 2>&1
    else
        echo "[$(date)] $srr already exists, skipping." >> "$LOG"
    fi
done

echo "[$(date)] All 28 SRA downloads finished successfully." >> "$LOG"
