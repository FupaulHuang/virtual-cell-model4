#!/usr/bin/env python3
"""Check a source-only checkout or staged Git contents without printing secrets."""
from pathlib import Path, PurePosixPath
import argparse
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
ROOT_FILES = {'README.md', 'PUBLICATION.md', 'THIRD_PARTY_NOTICES.md', 'PLAN.md',
              'EXPERIMENT_PLAN.md', 'SUBMISSION_PLAN.md', 'requirements.txt', '.gitignore'}
PATTERNS = {
    'credential': re.compile(r'\b(?:vcc_pat_|github_pat_|gh[pousr]_)[A-Za-z0-9_]{12,}'),
    'private_key': re.compile(r'-----BEGIN (?:[A-Z]+ )?PRIVATE KEY-----'),
    'machine_path': re.compile(r'/(?:home|mnt|Users)/[^\s]+'),
    'signed_url': re.compile(r'(?:X-Goog-Signature|X-Amz-Signature|upload_id|access_token)=[A-Za-z0-9%_-]{8,}', re.I),
    'email': re.compile(r'[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.(?:com|net|org|edu|gov|io|cn|hk|uk)\b', re.I),
    'bearer_value': re.compile(r'Bearer\s+[A-Za-z0-9._~-]{16,}'),
}


def permitted(name):
    p = PurePosixPath(name)
    if name in ROOT_FILES or name == '.githooks/pre-commit':
        return True
    if len(p.parts) == 2:
        folder, leaf = p.parts
        if folder in ('scripts', 'tests'):
            return p.suffix in ({'.py', '.sh'} if folder == 'scripts' else {'.py'})
        if folder == 'configs':
            return p.suffix == '.json'
        if folder in ('docs', 'reports', 'research'):
            return p.suffix == '.md' or name == 'reports/upstream_sources.json'
    if len(p.parts) == 3 and p.parts[:2] == ('external', 'atlasshift'):
        return p.suffix == '.py' or p.name in ('LICENSE', 'README.md', 'requirements.txt', 'sources.json')
    return False


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--staged', action='store_true', help='Inspect the complete index, not unstaged working files.')
    args = ap.parse_args()
    if args.staged:
        names = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode().split('\0')
        files = [(n, subprocess.check_output(['git', 'show', ':' + n], cwd=ROOT)) for n in names if n]
    else:
        files = [(p.relative_to(ROOT).as_posix(), p.read_bytes()) for p in ROOT.rglob('*')
                 if p.is_file() and '.git' not in p.relative_to(ROOT).parts
                 and '__pycache__' not in p.parts and '.pytest_cache' not in p.parts]
    failures = []
    for name, raw in files:
        if not permitted(name):
            failures.append((name, 'unapproved_file'))
            continue
        if len(raw) > 2 * 1024**2:
            failures.append((name, 'oversize_source_file'))
        try:
            text = raw.decode('utf-8')
        except UnicodeDecodeError:
            failures.append((name, 'binary_content'))
            continue
        for kind, pattern in PATTERNS.items():
            if pattern.search(text):
                failures.append((name, kind))
    if failures:
        for name, kind in failures:
            print(f'FAIL {name}: {kind}', file=sys.stderr)
        raise SystemExit(1)
    print(f'PASS: {len(files)} allowlisted text files; no configured sensitive patterns found.')


if __name__ == '__main__':
    main()
