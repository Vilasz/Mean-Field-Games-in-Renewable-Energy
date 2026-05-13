"""Quick inspection helper for notebooks 04 and 05."""
import sys
import nbformat


def show(path: str) -> None:
    nb = nbformat.read(path, as_version=4)
    print(f"=== {path}  ({len(nb.cells)} cells) ===")
    for i, c in enumerate(nb.cells):
        head = c.source.replace("\n", " / ")[:160]
        print(f"{i:>3}: [{c.cell_type:8}] {head}")


if __name__ == "__main__":
    for p in sys.argv[1:]:
        show(p)
