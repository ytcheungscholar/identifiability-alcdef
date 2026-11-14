#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Fill (or check) the submission placeholders in one place.

Usage
-----
    python fill_placeholders.py --check
    python fill_placeholders.py --email you@uni.edu --doi 10.5281/zenodo.1234567 \
                               --github your-username [--name "Yuen Tai Cheung"]

`--check` only reports what is still a placeholder and exits non-zero if any
remain, so it can be run right before building the final PDF.
"""

import argparse
import io
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.join(os.path.dirname(HERE), "paper_ac_20261114")
REPO = HERE


def targets():
    return [
        (os.path.join(PAPER, "main.tex"), "paper"),
        (os.path.join(PAPER, "cover_letter.md"), "paper"),
        (os.path.join(REPO, "CITATION.cff"), "repo"),
        (os.path.join(REPO, "README.md"), "repo"),
        (os.path.join(REPO, "LICENSE"), "repo"),
    ]


PLACEHOLDER_PATTERNS = [
    (r"your\.email@example\.com", "email"),
    (r"zenodo\.XXXXXXX", "zenodo DOI"),
    (r"\bXXXXX\b", "GitHub account"),
    (r"<your-username>", "GitHub account"),
]


def scan(verbose=True):
    pending = []
    for path, kind in targets():
        if not os.path.exists(path):
            continue
        text = io.open(path, encoding="utf-8-sig").read()
        for lineno, line in enumerate(text.splitlines(), 1):
            for pat, label in PLACEHOLDER_PATTERNS:
                if re.search(pat, line):
                    pending.append((path, lineno, label, line.strip()))
                    if verbose:
                        print("%-16s %-46s L%-4d %s"
                              % (label, os.path.basename(path), lineno, line.strip()[:78]))
    return pending


def apply(email=None, doi=None, github=None, name=None):
    for path, kind in targets():
        if not os.path.exists(path):
            continue
        text = io.open(path, encoding="utf-8-sig").read()
        orig = text
        if email:
            text = text.replace("your.email@example.com", email)
        if doi:
            text = text.replace("zenodo.XXXXXXX", doi.replace("https://doi.org/", ""))
        if github:
            text = text.replace("/XXXXX/", "/%s/" % github).replace("<your-username>", github)
        if name:
            text = text.replace("Yuen Tai Cheung", name)
        if text != orig:
            io.open(path, "w", encoding="utf-8", newline="\n").write(text)
            print("patched:", os.path.basename(path))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--email")
    ap.add_argument("--doi")
    ap.add_argument("--github")
    ap.add_argument("--name")
    a = ap.parse_args()

    if a.check:
        print("剩餘占位符：\n")
        pending = scan()
        if not pending:
            print("(none) — 可以編最終版了")
            return 0
        print("\n共 %d 處，請用 --email/--doi/--github 填寫或手動編輯" % len(pending))
        return 1

    apply(a.email, a.doi, a.github, a.name)
    print("\n現在狀態：")
    pending = scan()
    return 1 if pending else 0


if __name__ == "__main__":
    sys.exit(main())
