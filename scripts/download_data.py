import os
import urllib.request
import time

URL = "https://opendata.cern.ch/record/303/files/dimuon.csv"
OUT_PATH = os.path.join(os.path.dirname(__file__), "..", "data", "lhc", "dimuon.csv")

def download_file():
    print(f"Iniciando download de {URL}...")
    req = urllib.request.Request(
        URL,
        headers={"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36"}
    )
    start = time.time()
    with urllib.request.urlopen(req, timeout=60) as resp, open(OUT_PATH, "wb") as f:
        total = int(resp.headers.get("Content-Length", 0))
        downloaded = 0
        chunk_size = 1024 * 1024
        while True:
            chunk = resp.read(chunk_size)
            if not chunk:
                break
            f.write(chunk)
            downloaded += len(chunk)
            if total:
                print(f"Progresso: {downloaded / 1024 / 1024:.2f} MB / {total / 1024 / 1024:.2f} MB ({downloaded/total*100:.1f}%)")
    elapsed = time.time() - start
    print(f"Download concluido em {elapsed:.1f}s: {OUT_PATH} ({os.path.getsize(OUT_PATH):,} bytes)")

if __name__ == "__main__":
    download_file()
