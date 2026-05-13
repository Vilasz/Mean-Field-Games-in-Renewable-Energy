"""Executa um notebook in-place e retorna 0 se OK."""
from __future__ import annotations

import sys
import time
from pathlib import Path

import nbformat
from nbclient import NotebookClient


def run(path_str: str, kernel_name: str = "python310", timeout: int = 1800) -> int:
    p = Path(path_str)
    nb = nbformat.read(p, as_version=4)
    client = NotebookClient(
        nb,
        timeout=timeout,
        kernel_name=kernel_name,
        resources={"metadata": {"path": str(p.parent)}},
    )
    t0 = time.perf_counter()
    try:
        client.execute()
    except Exception as exc:
        print(f"[FAIL] {p.name}: {type(exc).__name__}: {exc}")
        nbformat.write(nb, p)
        return 1
    nbformat.write(nb, p)
    dt = time.perf_counter() - t0
    print(f"[OK] {p.name} executado em {dt:.1f}s")
    return 0


if __name__ == "__main__":
    rcs = [run(arg) for arg in sys.argv[1:]]
    sys.exit(max(rcs) if rcs else 0)
