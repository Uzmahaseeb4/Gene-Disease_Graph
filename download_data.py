#!/usr/bin/env python3
"""
Download raw files for the Gene-Disease KG project (no wget needed, stdlib only).

Run from your project folder:
    python download_data.py
Creates:
    data/raw/opentargets/{target,disease,association_overall_direct}/*.parquet
    data/raw/gwas/gwas-catalog-associations_ontology-annotated.tsv

If a URL returns 404, the release layout changed: open the base URL in a
browser, check the folder names, and edit OT_BASE / OT_DATASETS below.
"""
import os
import re
import sys
import time
import zipfile
import urllib.request

OT_BASE = "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/latest/output/"
OT_DATASETS = ["target", "disease", "association_overall_direct"]
GWAS_DIR = "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/"
RAW = os.path.join("data", "raw")


def remote_size(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return int(r.headers.get("Content-Length", -1))


def download(url, dest, retries=15, chunk=1024 * 1024):
    """Resumable download: writes to dest + '.part', renames when complete."""
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    total = remote_size(url)
    if os.path.exists(dest) and (total < 0 or os.path.getsize(dest) == total):
        print(f"  skip (complete): {dest}")
        return
    part = dest + ".part"
    for attempt in range(1, retries + 1):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        if total > 0 and have == total:
            break
        req = urllib.request.Request(url, headers={"Range": f"bytes={have}-"} if have else {})
        try:
            with urllib.request.urlopen(req, timeout=60) as r:
                if have and r.status != 206:      # server ignored Range -> restart
                    have = 0
                with open(part, "ab" if have else "wb") as f:
                    while True:
                        buf = r.read(chunk)
                        if not buf:
                            break
                        f.write(buf)
            have = os.path.getsize(part)
            if total < 0 or have == total:
                break
            print(f"  incomplete ({have}/{total}), retry {attempt}/{retries}")
        except Exception as e:
            print(f"  network error ({e!r}), retry {attempt}/{retries}")
            time.sleep(min(2 * attempt, 20))
    if total > 0 and os.path.getsize(part) != total:
        raise RuntimeError(f"Could not finish {url}; rerun the script to resume.")
    os.replace(part, dest)
    print(f"  done: {dest}")


def list_parquet(url):
    html = urllib.request.urlopen(url).read().decode("utf-8", "ignore")
    return sorted(set(re.findall(r'href="([^"?/][^"]*\.parquet)"', html)))


def get_datatype_dataset():
    """Associations split by evidence TYPE (genetic, drug, literature, ...). The folder name varies
    between Open Targets releases, so find it by looking for 'datatype' + 'direct' (not 'indirect')."""
    print("[Open Targets] association by datatype (direct)")
    html = urllib.request.urlopen(OT_BASE, timeout=60).read().decode("utf-8", "ignore")
    names = sorted(set(re.findall(r'href="([^"?/][^"]*)/"', html)))
    cands = [n for n in names if "datatype" in n.lower() and "direct" in n.lower()
             and "indirect" not in n.lower()]
    print("  association-related folders:", [n for n in names if "assoc" in n.lower()])
    if not cands:
        print("  No datatype/direct folder found. Pick one from the list above and add it to OT_DATASETS.")
        return
    name = cands[0]
    url = f"{OT_BASE}{name}/"
    files = list_parquet(url)
    print(f"  using: {name} ({len(files)} files)")
    for f in files:
        download(url + f, os.path.join(RAW, "opentargets", name, f))


def get_gwas():
    """Find the current 'all associations' file on the EBI FTP, download, unzip."""
    html = urllib.request.urlopen(GWAS_DIR, timeout=60).read().decode("utf-8", "ignore")
    names = sorted(set(re.findall(r'href="([^"?/][^"]*)"', html)))
    cands = [n for n in names
             if "associations" in n.lower() and n.lower().endswith((".zip", ".tsv"))
             and not re.search(r"split|year|summary|studies|ancestry", n, re.I)]
    # need MAPPED_TRAIT_URI -> prefer the ontology-annotated file
    cands.sort(key=lambda n: ("ontology-annotated" not in n, not n.endswith(".zip"), n))
    print("  candidates:", cands)
    if not cands:
        print("  No matching file. Files in folder:", names)
        print("  Download manually from https://www.ebi.ac.uk/gwas/docs/file-downloads")
        return
    pick = cands[0]
    print("  using:", pick)
    dest = os.path.join(RAW, "gwas", pick)
    download(GWAS_DIR + pick, dest)
    if pick.endswith(".zip"):
        with zipfile.ZipFile(dest) as z:
            for member in z.namelist():
                if member.lower().endswith(".tsv"):
                    z.extract(member, os.path.join(RAW, "gwas"))
                    print("  extracted:", member)
        os.remove(dest)


def main():
    for name in OT_DATASETS:
        url = f"{OT_BASE}{name}/"
        print(f"[Open Targets] {name}")
        try:
            files = list_parquet(url)
        except Exception as e:
            print(f"  ERROR listing {url}: {e}")
            continue
        if not files:
            print("  no parquet files found; check the URL in a browser")
        for f in files:
            download(url + f, os.path.join(RAW, "opentargets", name, f))

    get_datatype_dataset()

    print("[GWAS Catalog]")
    get_gwas()
    print("Done.")


if __name__ == "__main__":
    sys.exit(main())
