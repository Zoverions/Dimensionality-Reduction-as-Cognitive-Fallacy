import requests
import numpy as np
import time
import os

# --- CONFIGURATION ---
API_KEY = "9d963b728d2ffbc4778b7b3f5188e742"
HEADERS = {"api-key": API_KEY}
BASE_URL = "https://www.tng-project.org/api/TNG300-1/snapshots/"
SNAPSHOTS = {"z=0": 99, "z=1": 50}
LIMIT = 15000

def fetch_robust_catalog(snap_num, limit):
    print(f"--- Fetching Top {limit} Subhalos for Snapshot {snap_num} ---")
    
    # First, get all subhalo IDs and masses from the list endpoint
    ids_url = f"{BASE_URL}{snap_num}/subhalos/"
    params = {"limit": limit, "order_by": "-mass_log_msun", "fields": "id,mass_log_msun"}
    
    ids, masses_log = [], []
    next_url = ids_url
    while next_url and len(ids) < limit:
        try:
            r = requests.get(next_url, params=params if next_url == ids_url else {}, headers=HEADERS)
            if r.status_code == 429:
                time.sleep(int(r.headers.get("Retry-After", 5)))
                continue
            r.raise_for_status()
            data = r.json()
            for sub in data["results"]:
                ids.append(sub["id"])
                masses_log.append(sub["mass_log_msun"])
            print(f"  Fetched {len(ids)} / {limit} IDs...", end="\r")
            next_url = data.get("next")
        except Exception as e:
            print(f"\nAPI Error fetching IDs: {e}")
            return None, None

    print(f"\nFetched {len(ids)} IDs. Now fetching positions...")

    # Now, fetch positions for each ID individually
    positions = []
    for i, sub_id in enumerate(ids):
        pos_url = f"{BASE_URL}{snap_num}/subhalos/{sub_id}/"
        try:
            r = requests.get(pos_url, headers=HEADERS)
            if r.status_code == 429:
                time.sleep(int(r.headers.get("Retry-After", 5)))
                # Retry the same ID
                r = requests.get(pos_url, headers=HEADERS)
            r.raise_for_status()
            data = r.json()
            positions.append([data["pos_x"], data["pos_y"], data["pos_z"]])
            if (i + 1) % 100 == 0:
                print(f"  Fetched {i+1} / {len(ids)} positions...", end="\r")
        except Exception as e:
            print(f"\nAPI Error fetching position for ID {sub_id}: {e}")
            # Pad with zeros to maintain array shape
            positions.append([0, 0, 0])

    print(f"\nComplete. Retrieved {len(positions)} positions.")

    # Convert to numpy arrays
    pos_arr = np.array(positions, dtype=np.float64)
    mass_arr = 10**np.array(masses_log, dtype=np.float64) # Convert log mass to linear

    # Unit Conversion: ckpc/h -> cMpc/h (Comoving)
    pos_arr /= 1000.0

    return pos_arr, mass_arr

if __name__ == "__main__":
    for label, snap in SNAPSHOTS.items():
        pos, mass = fetch_robust_catalog(snap, LIMIT)
        if pos is not None and mass is not None:
            np.save(f"tng300_{label}_pos_{LIMIT}.npy", pos)
            np.save(f"tng300_{label}_mass_{LIMIT}.npy", mass)
            print(f"Saved {label}: Pos shape {pos.shape}, Mass shape {mass.shape}")
