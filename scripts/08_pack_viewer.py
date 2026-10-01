"""Copy the export into the viewer. edges.bin is shipped as base64 inside a .json file,
because download managers (e.g. IDM) intercept .bin / octet-stream requests."""
import base64, json, shutil
from pathlib import Path
src, dst = Path("data/export"), Path("viewer/public/data")
dst.mkdir(parents=True, exist_ok=True)
for f in src.iterdir():
    if f.name != "edges.bin":
        shutil.copy(f, dst / f.name)
raw = (src / "edges.bin").read_bytes()
(dst / "edges_b64.json").write_text(json.dumps(base64.b64encode(raw).decode("ascii")), encoding="utf-8")
(dst / "edges.bin").unlink(missing_ok=True)
print(f"edges: {len(raw):,} bytes -> edges_b64.json")
