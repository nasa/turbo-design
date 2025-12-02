#!/usr/bin/env python3
"""Scan the repository for Windows-incompatible file names and fix them.

This script recursively searches for files and directories containing
characters that are illegal on Windows (e.g. *, ?, <, >, \\) or control
characters, as well as names that end with spaces/dots. Offending
characters are replaced with '-' and names that would collide are made
unique by appending a numeric suffix.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Iterable, List, Sequence, Tuple

INVALID_CHARS = set('<>:"/\\|?*#')
CONTROL_MAX = 31
RESERVED_NAMES = {
    'CON',
    'PRN',
    'AUX',
    'NUL',
    *{f'COM{i}' for i in range(1, 10)},
    *{f'LPT{i}' for i in range(1, 10)},
}
DEFAULT_SKIPS = {'.git', '.hg', '.svn', '.tox'}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Replace Windows-incompatible characters in file and directory names.'
    )
    parser.add_argument(
        '--root',
        type=Path,
        default=Path.cwd(),
        help='Root directory to scan (defaults to current working directory).',
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show planned renames without applying them.',
    )
    parser.add_argument(
        '--include-hidden',
        action='store_true',
        help='Include hidden files/directories (those that start with .).',
    )
    parser.add_argument(
        '--skip',
        action='append',
        default=[],
        help='Additional directory names to skip during the scan.',
    )
    return parser.parse_args()


def should_skip(path: Path, include_hidden: bool, skip_names: Sequence[str]) -> bool:
    name = path.name
    if name in skip_names:
        return True
    if not include_hidden and name.startswith('.'):
        return True
    return False


def sanitize_component(name: str) -> Tuple[str, bool]:
    changed = False
    sanitized_chars: List[str] = []
    for ch in name:
        if ord(ch) <= CONTROL_MAX or ch in INVALID_CHARS:
            sanitized_chars.append('-')
            changed = True
        else:
            sanitized_chars.append(ch)
    sanitized = ''.join(sanitized_chars)

    trimmed = sanitized.rstrip(' .')
    if trimmed != sanitized:
        sanitized = trimmed
        changed = True

    if not sanitized:
        sanitized = '-'
        changed = True

    if sanitized.upper() in RESERVED_NAMES:
        sanitized = f'{sanitized}-'
        changed = True

    return sanitized, changed


def dedupe_name(parent: Path, candidate: str) -> Tuple[str, bool]:
    if not (parent / candidate).exists():
        return candidate, False

    stem, suffix = os.path.splitext(candidate)
    counter = 1
    while True:
        candidate_name = f'{stem}-{counter}{suffix}'
        target = parent / candidate_name
        if not target.exists():
            return candidate_name, True
        counter += 1


def collect_paths(root: Path, include_hidden: bool, skip_names: Sequence[str]) -> Tuple[List[Path], List[Path]]:
    file_paths: List[Path] = []
    dir_paths: List[Path] = []
    for current_root, dirs, files in os.walk(root):
        current_path = Path(current_root)
        filtered_dirs = []
        for directory in dirs:
            dir_path = current_path / directory
            if should_skip(dir_path, include_hidden, skip_names):
                continue
            filtered_dirs.append(directory)
            dir_paths.append(dir_path)
        dirs[:] = filtered_dirs

        for file_name in files:
            file_path = current_path / file_name
            if should_skip(file_path, include_hidden, skip_names):
                continue
            file_paths.append(file_path)

    dir_paths.sort(key=lambda p: len(p.relative_to(root).parts), reverse=True)
    return file_paths, dir_paths


def rename_path(path: Path, dry_run: bool) -> Tuple[bool, str]:
    parent = path.parent
    new_name, changed = sanitize_component(path.name)
    if not changed:
        return False, ''

    unique_name, deduped = dedupe_name(parent, new_name)
    target = parent / unique_name

    if dry_run:
        return True, f'{path} -> {target}'

    path.rename(target)
    action_note = ' (deduped)' if deduped else ''
    return True, f'{path} -> {target}{action_note}'


def process_paths(paths: Iterable[Path], dry_run: bool) -> List[str]:
    results: List[str] = []
    for path in paths:
        changed, message = rename_path(path, dry_run)
        if changed:
            results.append(message)
    return results


def main() -> None:
    args = parse_args()
    root = args.root.resolve()
    if not root.exists():
        raise SystemExit(f'Root path {root} does not exist')

    skip_names = set(DEFAULT_SKIPS)
    skip_names.update(args.skip)

    files, dirs = collect_paths(root, args.include_hidden, skip_names)
    file_changes = process_paths(files, args.dry_run)
    dir_changes = process_paths(dirs, args.dry_run)

    total_changes = len(file_changes) + len(dir_changes)
    if total_changes == 0:
        print('All file and directory names are already Windows-safe.')
        return

    print('\n'.join(file_changes + dir_changes))
    summary = 'Planned' if args.dry_run else 'Applied'
    print(f'{summary} {total_changes} rename(s).')


if __name__ == '__main__':
    main()
