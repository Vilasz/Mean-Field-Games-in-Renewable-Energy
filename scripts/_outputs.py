"""Print outputs of selected cells from executed notebooks."""
from __future__ import annotations

import json
import sys
from pathlib import Path

import nbformat


def show_text_outputs(path_str: str, only_idx: list[int] | None = None) -> None:
    nb = nbformat.read(Path(path_str), as_version=4)
    print(f"\n========================= {path_str} =========================")
    for i, c in enumerate(nb.cells):
        if c.cell_type != "code":
            continue
        if only_idx is not None and i not in only_idx:
            continue
        outputs = c.get("outputs", [])
        text = []
        for o in outputs:
            if o.output_type == "stream":
                text.append(o.text)
            elif o.output_type in ("execute_result", "display_data"):
                d = o.get("data", {})
                if "text/plain" in d:
                    text.append(d["text/plain"])
            elif o.output_type == "error":
                text.append("ERROR: " + " ".join(o.get("traceback", [])))
        if text:
            print(f"\n--- cell {i} ---")
            print("".join(text)[:4000])


if __name__ == "__main__":
    # NB4: capture audit (idx 9) and cota regression (idx 5), summary (idx 19)
    show_text_outputs("04_despacho_hidrotermico.ipynb")
    show_text_outputs("05_prob_diario_cmo_termico.ipynb")
