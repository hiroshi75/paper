#!/usr/bin/env python3
"""
Download GBIF occurrence data and WorldClim climate data for the
climate refugia thermal asymmetry analysis.

This script downloads:
1. GBIF occurrence records for 9 European tree/herb species
2. WorldClim 2.1 bioclimatic variables at 10 arc-minute resolution

Usage:
    python download_data.py [--output-dir OUTPUT_DIR] [--skip-worldclim]
"""

import os
import json
import time
import argparse
import zipfile
import requests
import csv

# Study region: Western Europe
LAT_MIN, LAT_MAX = 36.0, 55.0
LON_MIN, LON_MAX = -10.0, 15.0
MAX_RECORDS = 5000

# Species list with GBIF taxon keys
SPECIES = {
    "Quercus robur": 2878688,
    "Robinia pseudoacacia": 5352251,
    "Fagus sylvatica": 2882316,
    "Ailanthus altissima": 3190653,
    "Picea abies": 5284884,
    "Prunus serotina": 3021974,
    "Quercus ilex": 2878833,
    "Phytolacca americana": 3084850,
    "Impatiens glandulifera": 3189846,
}

WORLDCLIM_URL = "https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_10m_bio.zip"


def download_gbif_species(species_name, taxon_key, output_dir):
    """Download occurrence records for a single species from GBIF API."""
    safe_name = species_name.lower().replace(" ", "_")
    output_file = os.path.join(output_dir, f"gbif_{safe_name}.csv")

    if os.path.exists(output_file):
        print(f"  [SKIP] {output_file} already exists")
        return output_file

    print(f"  Downloading {species_name} (taxon key: {taxon_key})...")

    records = []
    offset = 0
    limit = 300

    while len(records) < MAX_RECORDS:
        params = {
            "taxonKey": taxon_key,
            "hasCoordinate": "true",
            "hasGeospatialIssue": "false",
            "decimalLatitude": f"{LAT_MIN},{LAT_MAX}",
            "decimalLongitude": f"{LON_MIN},{LON_MAX}",
            "limit": limit,
            "offset": offset,
        }
        url = "https://api.gbif.org/v1/occurrence/search"

        try:
            resp = requests.get(url, params=params, timeout=60)
            resp.raise_for_status()
            data = resp.json()
        except Exception as e:
            print(f"    Warning: API error at offset {offset}: {e}")
            break

        results = data.get("results", [])
        if not results:
            break

        for r in results:
            lat = r.get("decimalLatitude")
            lon = r.get("decimalLongitude")
            year = r.get("year")
            if lat is not None and lon is not None:
                records.append({
                    "species": species_name,
                    "decimalLatitude": lat,
                    "decimalLongitude": lon,
                    "year": year,
                    "gbifID": r.get("gbifID"),
                })

        offset += limit
        if not data.get("endOfRecords", True) is False:
            break

        time.sleep(0.5)  # Rate limiting

    # Write CSV
    if records:
        with open(output_file, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["species", "decimalLatitude",
                                                     "decimalLongitude", "year", "gbifID"])
            writer.writeheader()
            writer.writerows(records)
        print(f"    Saved {len(records)} records to {output_file}")
    else:
        print(f"    Warning: No records retrieved for {species_name}")

    return output_file


def download_worldclim(output_dir):
    """Download WorldClim 2.1 bioclimatic variables at 10 arc-min resolution."""
    wc_dir = os.path.join(output_dir, "worldclim")
    zip_path = os.path.join(output_dir, "wc2.1_10m_bio.zip")

    # Check if already extracted
    expected_file = os.path.join(wc_dir, "wc2.1_10m_bio_1.tif")
    if os.path.exists(expected_file):
        print(f"  [SKIP] WorldClim data already exists at {wc_dir}")
        return wc_dir

    os.makedirs(wc_dir, exist_ok=True)

    print(f"  Downloading WorldClim 2.1 (~90 MB)...")
    try:
        resp = requests.get(WORLDCLIM_URL, stream=True, timeout=300)
        resp.raise_for_status()
        with open(zip_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)
        print(f"    Downloaded to {zip_path}")
    except Exception as e:
        print(f"    Error downloading WorldClim: {e}")
        print(f"    Please download manually from: {WORLDCLIM_URL}")
        print(f"    Extract to: {wc_dir}")
        return wc_dir

    # Extract
    print(f"  Extracting...")
    try:
        with zipfile.ZipFile(zip_path, "r") as z:
            z.extractall(wc_dir)
        os.remove(zip_path)
        print(f"    Extracted to {wc_dir}")
    except Exception as e:
        print(f"    Error extracting: {e}")

    return wc_dir


def main():
    parser = argparse.ArgumentParser(description="Download data for refugia analysis")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--skip-worldclim", action="store_true",
                        help="Skip WorldClim download")
    args = parser.parse_args()

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("Data Download for Climate Refugia Thermal Asymmetry Analysis")
    print("=" * 60)

    # 1. GBIF data
    print("\n[1/2] Downloading GBIF occurrence data...")
    for species_name, taxon_key in SPECIES.items():
        download_gbif_species(species_name, taxon_key, output_dir)
        time.sleep(1)

    # 2. WorldClim data
    if not args.skip_worldclim:
        print("\n[2/2] Downloading WorldClim climate data...")
        download_worldclim(output_dir)
    else:
        print("\n[2/2] Skipping WorldClim download")

    print("\n" + "=" * 60)
    print("Download complete!")
    print(f"Data saved to: {os.path.abspath(output_dir)}")
    print("=" * 60)


if __name__ == "__main__":
    main()
