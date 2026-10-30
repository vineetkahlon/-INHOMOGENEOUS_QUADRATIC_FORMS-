"""Detect changes to delivered files using the supplied SHA-256 manifest."""
import hashlib,json
from pathlib import Path

def main():
    root=Path(__file__).resolve().parent
    manifest=json.loads((root/'manifest_sha256.json').read_text())
    for name,expected in manifest.items():
        path=root/name
        if not path.is_file():raise ValueError(f'Missing file: {name}')
        if hashlib.sha256(path.read_bytes()).hexdigest()!=expected:
            raise ValueError(f'Hash mismatch: {name}')
    print(json.dumps({'files_checked':len(manifest),'result':'All listed file hashes match.'},indent=2))
if __name__=='__main__':main()
