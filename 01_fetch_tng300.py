"""
01_fetch_tng300.py
Fetches the top 15,000 most massive subhalos from IllustrisTNG300-1 for z=1 and z=0.
Requires an API key from www.tng-project.org.
"""

import requests
import numpy as np
import os

# --- CONFIGURATION ---
API_KEY = 'YOUR_API_KEY_HERE' # Replace with your TNG API key
BASE_URL = 'http://www.tng-project.org/api/TNG300-1/'
HEADERS = {"api-key": API_KEY}
N_TARGET = 15000

# Snapshots for z=1 and z=0 in TNG300-1
SNAPSHOTS = {'z=1': 50, 'z=0': 99}

def get_subhalo_catalog(snapshot, n_target):
    """Fetches mass and position for the most massive subhalos."""
    print(f"Fetching data for Snapshot {snapshot}...")

    # We request specific fields to minimize payload size
    search_query = f"?limit={n_target}&order_by=-SubhaloMass&subhalos_fields=SubhaloMass,SubhaloPos"
    url = f"{BASE_URL}snapshots/{snapshot}/subhalos/{search_query}"

    response = requests.get(url, headers=HEADERS)
    response.raise_for_status()
    data = response.json()

    masses = []
    positions = []

    # The API paginates, but a large limit usually covers the top N.
    # We iterate through the results to build our arrays.
    for subhalo in data['results']:
        masses.append(subhalo['SubhaloMass'])
        positions.append(subhalo['SubhaloPos'])

    return np.array(positions), np.array(masses)

if __name__ == "__main__":
    os.makedirs("data", exist_ok=True)

    for label, snap in SNAPSHOTS.items():
        try:
            pos, mass = get_subhalo_catalog(snap, N_TARGET)

            # Save to numpy arrays for the analysis script
            np.save(f"data/tng300_{label}_pos.npy", pos)
            np.save(f"data/tng300_{label}_mass.npy", mass)

            print(f"Successfully saved {len(mass)} halos for {label}.")

        except requests.exceptions.HTTPError as e:
            print(f"HTTP Error failed to fetch {label}: {e}")
            print("Did you insert your API key?")