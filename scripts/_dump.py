"""Dump full source of every cell in a notebook to stdout."""
import sys
import nbformat


def dump(path: str) -> None:
    nb = nbformat.read(path, as_version=4)
    print(f"\n========================= {path} =========================")
    for i, c in enumerate(nb.cells):
        print(f"\n--- cell {i:>3} [{c.cell_type}] ---")
        print(c.source)


if __name__ == "__main__":
    for p in sys.argv[1:]:
        dump(p)
