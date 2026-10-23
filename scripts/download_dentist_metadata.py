"""
PhyloGap - Script to download DENTIST benchmark gap coordinate metadata and ground truth reports.
Downloads Drosophila melanogaster and Calypte anna (Anna's hummingbird) benchmark definitions.
"""

import os
import sys
import urllib.request
import time

BASE_URL = "https://bds.mpi-cbg.de/hillerlab/DENTIST/data"
TARGET_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "benchmarks", "dentist_pilot")

DOWNLOAD_TARGETS = [
    # Drosophila melanogaster (Pilot model organism)
    {
        "organism": "d_melanogaster",
        "file": "assembly-reference.mapped-contigs.bed",
        "desc": "D. melanogaster mapped contig coordinates (test gaps are inter-contig regions)"
    },
    {
        "organism": "d_melanogaster",
        "file": "closable-gaps.assembly-reference.reads-simulated-pb.json",
        "desc": "D. melanogaster closable gaps report (simulated PacBio reads)"
    },
    # Calypte anna (Anna's hummingbird - Avian focal organism)
    {
        "organism": "c_anna",
        "file": "assembly-reference.mapped-contigs.bed",
        "desc": "C. anna (hummingbird) mapped contig coordinates"
    },
]

def download_file(url, target_path, desc):
    print(f"[*] Downloading: {desc}")
    print(f"    Source: {url}")
    print(f"    Target: {target_path}")
    
    start_time = time.time()
    try:
        def progress(count, block_size, total_size):
            percent = int(count * block_size * 100 / (total_size if total_size > 0 else 1))
            sys.stdout.write(f"\r    Progress: {percent}% ({count * block_size / 1024:.1f} KB)")
            sys.stdout.flush()

        urllib.request.urlretrieve(url, target_path, reporthook=progress)
        elapsed = time.time() - start_time
        size_kb = os.path.getsize(target_path) / 1024
        print(f"\n    [+] Done: {size_kb:.2f} KB in {elapsed:.2f}s")
        return True
    except Exception as e:
        print(f"\n    [-] Error downloading {url}: {e}")
        return False

def main():
    os.makedirs(TARGET_DIR, exist_ok=True)
    print(f"=== PhyloGap: Downloading Benchmark Metadata ===")
    print(f"Directory: {TARGET_DIR}\n")
    
    success_count = 0
    for item in DOWNLOAD_TARGETS:
        organism_dir = os.path.join(TARGET_DIR, item["organism"])
        os.makedirs(organism_dir, exist_ok=True)
        url = f"{BASE_URL}/{item['organism']}/{item['file']}"
        target_path = os.path.join(organism_dir, item["file"])
        
        if os.path.exists(target_path) and os.path.getsize(target_path) > 0:
            print(f"[i] Already exists: {item['file']} ({os.path.getsize(target_path)/1024:.2f} KB)")
            success_count += 1
            continue
            
        if download_file(url, target_path, item["desc"]):
            success_count += 1
            
    print(f"\n=== Download Summary: {success_count}/{len(DOWNLOAD_TARGETS)} files ready. ===")

if __name__ == "__main__":
    main()
