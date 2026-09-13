#!/usr/bin/env python3
"""Write hashes.txt: sha256 of every tracked source/doc file (sorted).

Excludes VCS, caches, logs, and hashes.txt itself. Run after editing sources;
CI and reviewers can re-run to confirm the tree matches.
"""
import hashlib, os

HERE = os.path.dirname(os.path.abspath(__file__))
SKIP_DIRS = {".git", "__pycache__", ".pytest_cache", ".github"}  # CI config excluded
SKIP_NAMES = {"hashes.txt"}
INCLUDE_EXT = {".py", ".md", ".cff", ".txt", ".yml", ".yaml"}


def main():
    rows = []
    for root, dirs, files in os.walk(HERE):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS]
        for name in files:
            if name in SKIP_NAMES:
                continue
            ext = os.path.splitext(name)[1]
            if ext not in INCLUDE_EXT:
                continue
            path = os.path.join(root, name)
            rel = os.path.relpath(path, HERE)
            with open(path, "rb") as fh:
                digest = hashlib.sha256(fh.read()).hexdigest()
            rows.append(f"{digest}  {rel}")
    rows.sort(key=lambda r: r.split("  ", 1)[1])
    with open(os.path.join(HERE, "hashes.txt"), "w") as fh:
        fh.write("\n".join(rows) + "\n")
    print(f"wrote hashes.txt with {len(rows)} entries")


if __name__ == "__main__":
    main()
