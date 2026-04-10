#!/usr/bin/env python3
# SPDX-FileCopyrightText: Copyright (c) 2026 NVIDIA CORPORATION & AFFILIATES. All rights reserved.
# SPDX-License-Identifier: Apache-2.0
"""
CI-only patches for Jupyter notebooks (papermill + nested Docker on self-hosted GPU runners).

Keeps committed notebooks aligned with upstream; run with CI_PATCH_NOTEBOOKS=1 before papermill.
Requires NOTEBOOK_DATA_DIR and HOST_DATA_DIR in the environment (set by GitHub Actions).
"""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path


def _source_to_text(source: list | str) -> str:
    if isinstance(source, str):
        return source
    return "".join(source)


def _text_to_source(text: str) -> list[str]:
    if not text:
        return []
    return text.splitlines(keepends=True)


DATA_DIR_BLOCK_OLD_VARIANTS = (
    (
        'import os\n\n'
        '# Set the path to the data directory \n'
        'os.environ[\'DATA_DIR\'] = "/ephemeral"\n'
        'print("✅ HOST_PATH exported to shell")'
    ),
    (
        'import os\n\n'
        '# Set the path to the data directory \n'
        'os.environ[\'DATA_DIR\'] = "/ephemeral"\n'
        'print("✅ DATA_DIR exported to shell")'
    ),
)

DATA_DIR_BLOCK_NEW = (
    'import os\n\n'
    'if os.environ.get("NOTEBOOK_DATA_DIR"):\n'
    '    os.environ["DATA_DIR"] = os.environ["NOTEBOOK_DATA_DIR"]\n'
    'else:\n'
    '    os.environ["DATA_DIR"] = "/ephemeral"\n'
    'print("✅ DATA_DIR=", os.environ["DATA_DIR"])'
)

DOWNLOAD_SCRIPTS = (
    "download_germline.sh",
    "download_pangenome.sh",
    "download_vep.sh",
)

QA_CELL = {
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [
        "%%sh\n",
        "\n",
        "# CI: marker substring for HTML-based QA (see blueprint-github-test-image).\n",
        'echo "NIST7035_TAAGGCGA_L001_R1_001.vcf"\n',
        'ls -lh "$DATA_DIR/out/NIST7035_TAAGGCGA_L001_R1_001.vcf"\n',
    ],
}


def patch_download_cell(text: str) -> str:
    for script in DOWNLOAD_SCRIPTS:
        old = f"sudo bash ../scripts/{script} $DATA_DIR"
        if old not in text:
            continue
        new = (
            'if [ "$(id -u)" -eq 0 ]; then\n'
            f'  bash ../scripts/{script} "$DATA_DIR"\n'
            "else\n"
            f'  sudo bash ../scripts/{script} "$DATA_DIR"\n'
            "fi\n"
        )
        return text.replace(old, new, 1)
    return text


def patch_docker_volumes(text: str) -> str:
    text = text.replace(
        "--volume $DATA_DIR/ref:/workdir",
        "--volume ${HOST_DATA_DIR:-$DATA_DIR}/ref:/workdir",
    )
    text = text.replace(
        "--volume $DATA_DIR:/workdir",
        "--volume ${HOST_DATA_DIR:-$DATA_DIR}:/workdir",
    )
    return text


def patch_notebook(path: Path) -> bool:
    data = json.loads(path.read_text(encoding="utf-8"))
    cells = data.get("cells", [])
    changed = False

    for cell in cells:
        if cell.get("cell_type") != "code":
            continue
        text = _source_to_text(cell.get("source", []))
        orig = text

        for old in DATA_DIR_BLOCK_OLD_VARIANTS:
            if old in text:
                text = text.replace(old, DATA_DIR_BLOCK_NEW, 1)
                break
        text = patch_download_cell(text)
        text = patch_docker_volumes(text)

        if text != orig:
            cell["source"] = _text_to_source(text)
            changed = True

    # Germline: insert QA cell after deepvariant (upstream has no marker in bcftools output).
    if path.name == "germline_wes.ipynb":
        already = any(
            'echo "NIST7035_TAAGGCGA_L001_R1_001.vcf"' in _source_to_text(c.get("source", []))
            for c in cells
            if c.get("cell_type") == "code"
        )
        if not already:
            for i, cell in enumerate(cells):
                if cell.get("cell_type") != "code":
                    continue
                t = _source_to_text(cell.get("source", []))
                if "pbrun deepvariant" in t and "--use-wes-model" in t:
                    cells.insert(i + 1, json.loads(json.dumps(QA_CELL)))
                    changed = True
                    break

    if changed:
        path.write_text(json.dumps(data, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    return changed


def main() -> int:
    if os.environ.get("CI_PATCH_NOTEBOOKS") != "1":
        return 0
    if not os.environ.get("NOTEBOOK_DATA_DIR") or not os.environ.get("HOST_DATA_DIR"):
        print(
            "ci_patch_notebooks: NOTEBOOK_DATA_DIR and HOST_DATA_DIR must be set",
            file=sys.stderr,
        )
        return 1
    if len(sys.argv) < 2:
        print("usage: ci_patch_notebooks.py <notebook.ipynb>", file=sys.stderr)
        return 1
    path = Path(sys.argv[1])
    if not path.is_file():
        print(f"ci_patch_notebooks: not found: {path}", file=sys.stderr)
        return 1
    patch_notebook(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
