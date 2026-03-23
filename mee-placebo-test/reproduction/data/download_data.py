#!/usr/bin/env python3
"""
Download data for the placebo test framework case studies.

This script downloads GBIF occurrence data for all six case studies:
1. Pollinator phenology (Apidae, Europe)
2. Dawn chorus (songbird species, Xeno-canto)
3. Plant traits (European tree species)
4. Nightlight-insect communities (Lepidoptera + Odonata)
5. Bird richness (Passeriformes + Strigiformes)
6. Forest fragmentation (tropical/temperate birds)

Usage:
    python download_data.py [--output-dir OUTPUT_DIR] [--case CASE_NUMBER]
"""

import os
import json
import time
import argparse
import csv
import requests

# GBIF API base URL
GBIF_API = "https://api.gbif.org/v1/occurrence/search"

# Case study species/taxa configurations
CASE_STUDIES = {
    1: {
        "name": "Pollinator phenology",
        "taxa": {
            "Apidae": {"familyKey": 4334, "continent": "EUROPE"},
        },
        "max_records": 50000,
    },
    2: {
        "name": "Dawn chorus (Xeno-canto)",
        "taxa": {
            "Turdus_merula": {"taxonKey": 9515886},
            "Erithacus_rubecula": {"taxonKey": 9517753},
            "Parus_major": {"taxonKey": 9596474},
            "Sylvia_atricapilla": {"taxonKey": 2493091},
            "Fringilla_coelebs": {"taxonKey": 9598506},
        },
        "max_records": 5000,
        "note": "Xeno-canto audio data must be downloaded separately from https://www.xeno-canto.org",
    },
    3: {
        "name": "Plant traits and drought",
        "taxa": {
            "European_trees": {"higherTaxonKey": 6, "continent": "EUROPE"},
        },
        "max_records": 10000,
        "note": "WOODIV trait data must be obtained separately from the WOODIV database",
    },
    4: {
        "name": "Nightlight and insect communities",
        "taxa": {
            "Lepidoptera": {"orderKey": 797, "continent": "EUROPE"},
            "Odonata": {"orderKey": 789, "continent": "EUROPE"},
        },
        "max_records": 20000,
    },
    5: {
        "name": "Bird richness and ALAN",
        "taxa": {
            "Passeriformes": {"orderKey": 729, "continent": "EUROPE"},
            "Strigiformes": {"orderKey": 734, "continent": "EUROPE"},
            "Caprimulgiformes": {"orderKey": 735, "continent": "EUROPE"},
        },
        "max_records": 8000,
    },
    6: {
        "name": "Forest fragmentation",
        "taxa": {
            "Aves_tropical": {"classKey": 212, "decimalLatitude": "-23.5,23.5"},
            "Poales": {"orderKey": 1369, "continent": "EUROPE"},
        },
        "max_records": 10000,
    },
}


def download_gbif_taxon(taxon_name, params, max_records, output_dir):
    """Download GBIF records for a given taxon."""
    safe_name = taxon_name.lower().replace(" ", "_")
    output_file = os.path.join(output_dir, f"gbif_{safe_name}.csv")

    if os.path.exists(output_file):
        print(f"  [SKIP] {output_file} already exists")
        return output_file

    print(f"  Downloading {taxon_name}...")

    records = []
    offset = 0
    limit = 300

    query_params = {
        "hasCoordinate": "true",
        "hasGeospatialIssue": "false",
        "limit": limit,
    }
    # Add taxon-specific params
    for k, v in params.items():
        if k not in ("continent",):
            query_params[k] = v
    if "continent" in params:
        query_params["continent"] = params["continent"]

    while len(records) < max_records:
        query_params["offset"] = offset

        try:
            resp = requests.get(GBIF_API, params=query_params, timeout=60)
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
            if lat is not None and lon is not None:
                records.append({
                    "species": r.get("species", taxon_name),
                    "decimalLatitude": lat,
                    "decimalLongitude": lon,
                    "year": r.get("year"),
                    "month": r.get("month"),
                    "day": r.get("day"),
                    "eventDate": r.get("eventDate"),
                    "gbifID": r.get("gbifID"),
                })

        offset += limit
        if data.get("endOfRecords", False):
            break

        time.sleep(0.5)

    if records:
        with open(output_file, "w", newline="") as f:
            fieldnames = ["species", "decimalLatitude", "decimalLongitude",
                          "year", "month", "day", "eventDate", "gbifID"]
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(records[:max_records])
        print(f"    Saved {min(len(records), max_records)} records to {output_file}")
    else:
        print(f"    Warning: No records retrieved for {taxon_name}")

    return output_file


def main():
    parser = argparse.ArgumentParser(description="Download data for placebo test case studies")
    parser.add_argument("--output-dir", default=".", help="Output directory")
    parser.add_argument("--case", type=int, default=0,
                        help="Case study number (1-6), or 0 for all")
    args = parser.parse_args()

    output_dir = args.output_dir
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("Data Download for Placebo Test Framework Case Studies")
    print("=" * 60)

    cases = [args.case] if args.case > 0 else range(1, 7)

    for case_num in cases:
        case = CASE_STUDIES[case_num]
        print(f"\n--- Case {case_num}: {case['name']} ---")
        if "note" in case:
            print(f"  Note: {case['note']}")

        for taxon_name, params in case["taxa"].items():
            download_gbif_taxon(taxon_name, params, case["max_records"], output_dir)
            time.sleep(1)

    print("\n" + "=" * 60)
    print("Download complete!")
    print(f"\nAdditional data sources (manual download required):")
    print("  - Xeno-canto recordings: https://www.xeno-canto.org")
    print("  - VIIRS nightlight: https://eogdata.mines.edu/products/vnl/")
    print("  - WorldClim: https://www.worldclim.org")
    print("  - MODIS NDVI: https://lpdaac.usgs.gov/products/mod13a3v061/")
    print("  - Hansen Forest Change: https://glad.earthengine.app/view/global-forest-change")
    print("=" * 60)


if __name__ == "__main__":
    main()
