#!/usr/bin/env python3
"""
push_kaggle_phylogap.py - Command-line Dispatcher & Monitor for Kaggle Multi-Account Jobs.

Allows dispatching PhyloGap experiment jobs to any of the 3 configured Kaggle accounts
(thung192, TThanh13, RaH1111) and checking kernel execution status / GPU quotas.

Usage:
    python scripts/push_kaggle_phylogap.py --list
    python scripts/push_kaggle_phylogap.py --status
    python scripts/push_kaggle_phylogap.py --push benchmark-sota-eval
    python scripts/push_kaggle_phylogap.py --push-all
    python scripts/push_kaggle_phylogap.py --download benchmark-sota-eval

Author: PhyloGap Research Team
Date: 2026-10-23
"""

import sys
import argparse
from pathlib import Path

# Add kaggle_runner_ui to sys.path
UI_DIR = Path(__file__).resolve().parents[1] / "kaggle_runner_ui"
if str(UI_DIR) not in sys.path:
    sys.path.insert(0, str(UI_DIR))

try:
    from runner.project_config import load_accounts, load_projects, accounts_by_id
    from runner.job_builder import prepare_project_jobs, push_job, status_job, download_job_output
    from runner.kaggle_api import get_account_gpu_status, validate_account
except ImportError as e:
    print(f"[ERROR] Could not load kaggle_runner_ui modules: {e}", file=sys.stderr)
    sys.exit(1)


def cmd_list(accounts, project):
    print("\n" + "="*80)
    print("PHYLOGAP KAGGLE JOBS & MULTI-ACCOUNT MATRIX")
    print("="*80)
    for job in project.get("jobs", []):
        acc = accounts.get(job["account_id"], {})
        user = acc.get("username", "Unknown")
        print(f"• Job Name:    {job['name']}")
        print(f"  Account:     {job['account_id']} (@{user})")
        print(f"  Description: {job.get('description', 'N/A')}")
        print(f"  Commands:    {'; '.join(job.get('commands', {}).get('run', []))}")
        print("-" * 80)


def cmd_status(accounts, project):
    print("\n" + "="*80)
    print("KAGGLE GPU QUOTA & KERNEL STATUS")
    print("="*80)
    acc_map = accounts_by_id()
    for acc_id, acc in acc_map.items():
        val_res = validate_account(acc)
        is_val = val_res.get("ok", False)
        msg = val_res.get("message", "")
        status_info = get_account_gpu_status(acc)
        print(f"Account: {acc_id} (@{acc.get('username')})")
        print(f"  Token Valid: {is_val} ({msg})")
        print(f"  Active GPU Slots: {status_info.get('active_gpu_count', 0)} / 2")
        print(f"  Running Kernels: {status_info.get('running_kernels', [])}")
        print("-" * 80)

    print("\n[JOB EXECUTION STATUS]")
    prepared = prepare_project_jobs(project, acc_map)
    for p in prepared:
        st = status_job(p["account"], p["job"], force_refresh=True)
        status_text = st.stdout.strip() if st.stdout else (st.stderr.strip() if st.stderr else "unknown")
        print(f"• {p['slug']} ({p['kernel_id']}): {status_text}")


def cmd_push(project, accounts, target_job_name=None):
    acc_map = accounts_by_id()
    prepared = prepare_project_jobs(project, acc_map)
    for p in prepared:
        if target_job_name and p["job"]["name"] != target_job_name and p["slug"] != target_job_name:
            continue
        print(f"\n[PUSH] Dispatching {p['slug']} to Kaggle as {p['kernel_id']}...")
        res = push_job(p["account"], p["job_dir"])
        print(f"       Result: {res}")


def cmd_download(project, accounts, target_job_name):
    acc_map = accounts_by_id()
    for job in project.get("jobs", []):
        if job["name"] == target_job_name or f"phylogap-{job['name']}" == target_job_name:
            acc = acc_map.get(job["account_id"])
            print(f"[DOWNLOAD] Downloading output for {job['name']} (@{acc['username']})...")
            out_dir, res = download_job_output(project["name"], acc, job)
            print(f"[SUCCESS] Downloaded to {out_dir}. Result: {res}")
            return
    print(f"[ERROR] Job '{target_job_name}' not found in project.")


def main():
    parser = argparse.ArgumentParser(description="PhyloGap Kaggle Runner CLI")
    parser.add_argument("--list", action="store_true", help="List all configured jobs and accounts")
    parser.add_argument("--status", action="store_true", help="Check Kaggle quota and kernel status")
    parser.add_argument("--push", type=str, help="Push specific job by name")
    parser.add_argument("--push-all", action="store_true", help="Push all configured jobs in parallel")
    parser.add_argument("--download", type=str, help="Download output artifacts of a completed job")
    args = parser.parse_args()

    acc_dict = accounts_by_id()
    projs = load_projects().get("projects", [])
    if not projs:
        print("[ERROR] No projects found in data/projects.yaml", file=sys.stderr)
        sys.exit(1)
    project = projs[0]

    if args.list or len(sys.argv) == 1:
        cmd_list(acc_dict, project)
    elif args.status:
        cmd_status(acc_dict, project)
    elif args.push:
        cmd_push(project, acc_dict, target_job_name=args.push)
    elif args.push_all:
        cmd_push(project, acc_dict)
    elif args.download:
        cmd_download(project, acc_dict, args.download)


if __name__ == "__main__":
    main()
