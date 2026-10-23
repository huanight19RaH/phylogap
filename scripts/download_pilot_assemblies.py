"""
PhyloGap - Script to download Drosophila melanogaster reference and test assemblies.
Target: data/benchmarks/dentist_pilot/d_melanogaster/
- assembly-reference.fasta.gz (~40.2 MB)
- assembly-test.fasta.gz (~41.2 MB)
"""

import os
import sys
import urllib.request
import time

BASE_URL = "https://bds.mpi-cbg.de/hillerlab/DENTIST/data/d_melanogaster"
TARGET_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "data", "benchmarks", "dentist_pilot", "d_melanogaster")

FILES_TO_DOWNLOAD = [
    {
        "file": "assembly-reference.fasta.gz",
        "desc": "Ground-truth reference assembly (FASTA gzipped, ~40.2 MB)"
    },
    {
        "file": "assembly-test.fasta.gz",
        "desc": "Test assembly with realistic copied gaps (FASTA gzipped, ~41.2 MB)"
    }
]

def download_file(url, target_path, desc):
    print(f"[*] Downloading: {desc}")
    print(f"    Source: {url}")
    print(f"    Target: {target_path}")
    
    start_time = time.time()
    try:
        def progress(count, block_size, total_size):
            if total_size > 0:
                percent = int(count * block_size * 100 / total_size)
                mb_downloaded = (count * block_size) / (1024 * 1024)
                total_mb = total_size / (1024 * 1024)
                sys.stdout.write(f"\r    Progress: {percent}% ({mb_downloaded:.1f}/{total_mb:.1f} MB)")
            else:
                kb_downloaded = (count * block_size) / 1024
                sys.stdout.write(f"\r    Progress: {kb_downloaded:.1f} KB downloaded")
            sys.stdout.flush()

        urllib.request.urlretrieve(url, target_path, reporthook=progress)
        elapsed = time.time() - start_time
        size_mb = os.path.getsize(target_path) / (1024 * 1024)
        print(f"\n    [+] Completed: {size_mb:.2f} MB in {elapsed:.2f}s ({size_mb/elapsed:.2f} MB/s)")
        return True
    except Exception as e:
        print(f"\n    [-] Error downloading {url}: {e}")
        return False

def main():
    os.makedirs(TARGET_DIR, exist_ok=True)
    print("=== PhyloGap: Downloading Pilot FASTA Assemblies ===")
    print(f"Target Directory: {TARGET_DIR}\n")
    
    success = 0
    for item in FILES_TO_DOWNLOAD:
        target_path = os.path.join(TARGET_DIR, item["file"])
        if os.path.exists(target_path) and os.path.getsize(target_path) > 1000000:
            print(f"[i] Already exists: {item['file']} ({os.path.getsize(target_path)/(1024*1024):.2f} MB)")
            success += 1
            continue
            
        url = f"{BASE_URL}/{item['file']}"
        if download_file(url, target_path, item["desc"]):
            success += 1
            
    print(f"\n=== Download Finished: {success}/{len(FILES_TO_DOWNLOAD)} assemblies ready. ===")

if __name__ == "__main__":
    main()
