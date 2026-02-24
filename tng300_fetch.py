"""
Phase 3: Definitive Data Fetch — TNG300-1 (z=0 vs z=1)
=======================================================
Fetches the top N most massive subhalos from TNG300-1 via the public REST API.

Adapted for the real TNG API structure:
  - List endpoint returns: id, mass_log_msun, sfr, url
  - Detail endpoint returns: pos_x, pos_y, pos_z, mass, ...
  - Positions in ckpc/h → converted to cMpc/h (comoving)

Methodology:
  1. Fetch top N subhalo IDs from list endpoint (sorted by -mass_log_msun)
  2. Fetch pos_x/y/z and mass from detail endpoints in parallel
  3. Save as .npy files for the analysis pipeline
"""

import os
import requests
import numpy as np
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# --- CONFIGURATION ---
API_KEY   = "9d963b728d2ffbc4778b7b3f5188e742"
HEADERS   = {"api-key": API_KEY}
BASE_URL  = "https://www.tng-project.org/api/TNG300-1/snapshots/"
SNAPSHOTS = {'z=0': 99, 'z=1': 50}
LIMIT     = 1000   # Sandbox-feasible; scale to 15000 on a cluster
MAX_WORKERS = 10
OUT_DIR   = '/home/ubuntu/Dimensionality-Reduction-as-Cognitive-Fallacy'

def fetch_top_ids(snap_num, limit):
    """Fetch the top `limit` subhalo IDs sorted by descending mass."""
    ids = []
    url = f"{BASE_URL}{snap_num}/subhalos/"
    params = {'limit': min(limit, 500), 'order_by': '-mass_log_msun'}
    
    while url and len(ids) < limit:
        try:
            r = requests.get(url, params=params, headers=HEADERS, timeout=60)
            if r.status_code == 429:
                print("  Rate limit hit. Sleeping 5s...")
                time.sleep(5)
                continue
            r.raise_for_status()
            data = r.json()
            for s in data['results']:
                ids.append(s['id'])
                if len(ids) >= limit:
                    break
            url = data.get('next')
            params = {}  # pagination URL includes params
        except Exception as e:
            print(f"  API Error during ID fetch: {e}")
            time.sleep(3)
            continue
    
    return ids[:limit]

def fetch_subhalo_detail(snap_num, subhalo_id):
    """Fetch pos_x/y/z and mass for a single subhalo from the detail endpoint."""
    url = f"{BASE_URL}{snap_num}/subhalos/{subhalo_id}/"
    for attempt in range(3):
        try:
            r = requests.get(url, headers=HEADERS, timeout=30)
            if r.status_code == 429:
                time.sleep(5)
                continue
            r.raise_for_status()
            d = r.json()
            # Positions in ckpc/h → cMpc/h
            pos = np.array([d['pos_x'], d['pos_y'], d['pos_z']]) / 1000.0
            mass = d['mass']  # 1e10 M_sun/h
            return pos, mass
        except Exception:
            time.sleep(1.5 * (attempt + 1))
    return None, None

def fetch_snapshot_data(snap_num, limit):
    """Fetch top `limit` subhalos: IDs from list, then positions from detail."""
    cache_file = os.path.join(OUT_DIR, f'tng300_snap{snap_num}_N{limit}.npz')
    if os.path.exists(cache_file):
        print(f"  Loading cached data from {cache_file}")
        data = np.load(cache_file)
        return data['pos'], data['mass']
    
    print(f"  Fetching top {limit} subhalo IDs for snapshot {snap_num}...")
    ids = fetch_top_ids(snap_num, limit)
    actual_ids = len(ids)
    print(f"  Got {actual_ids} IDs. Fetching positions ({MAX_WORKERS} workers)...")
    
    if actual_ids < limit:
        print(f"  WARNING: Only {actual_ids} subhalos available (target: {limit})")
    
    pos_list, mass_list = [None] * actual_ids, [None] * actual_ids
    with ThreadPoolExecutor(max_workers=MAX_WORKERS) as ex:
        future_to_idx = {
            ex.submit(fetch_subhalo_detail, snap_num, sid): i
            for i, sid in enumerate(ids)
        }
        done = 0
        for fut in as_completed(future_to_idx):
            i = future_to_idx[fut]
            pos, mass = fut.result()
            if pos is not None:
                pos_list[i]  = pos
                mass_list[i] = mass
            done += 1
            if done % 100 == 0:
                print(f"    {done}/{actual_ids} fetched...")
    
    # Filter out failed fetches
    valid = [(p, m) for p, m in zip(pos_list, mass_list) if p is not None]
    print(f"  Successfully fetched {len(valid)}/{actual_ids} subhalos.")
    
    pos_arr  = np.array([v[0] for v in valid])
    mass_arr = np.array([v[1] for v in valid])
    
    # Shape assertion
    assert pos_arr.ndim == 2 and pos_arr.shape[1] == 3, \
        f"Shape error: expected (N, 3), got {pos_arr.shape}"
    
    # Cache for re-runs
    np.savez(cache_file, pos=pos_arr, mass=mass_arr)
    print(f"  Cached to {cache_file}")
    
    return pos_arr, mass_arr

if __name__ == "__main__":
    for label, snap in SNAPSHOTS.items():
        print(f"\n=== {label} (Snapshot {snap}) ===")
        pos, mass = fetch_snapshot_data(snap, LIMIT)
        
        # Also save as .npy for compatibility with the analysis script
        np.save(os.path.join(OUT_DIR, f"tng300_{label}_pos.npy"), pos)
        np.save(os.path.join(OUT_DIR, f"tng300_{label}_mass.npy"), mass)
        
        print(f"  Saved {label}: Pos shape {pos.shape}, Mass range [{mass.min():.1f}, {mass.max():.1f}]")
    
    print("\n=== Data fetch complete ===")
