"""Inspect public GWPS HDF5 using an existing prefix plus explicit curl range files.

No network calls occur in Python. Missing byte ranges are emitted for direct curl.
"""
from pathlib import Path
import argparse, io, json
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / 'data/raw/gwps_subset'
PREFIX = ROOT / 'data/raw/K562_gwps_raw_singlecell_01.h5ad'
SIZE = 65830941948


class MissingRange(Exception):
    pass


class CachedFile(io.RawIOBase):
    def __init__(self):
        self.pos = 0
        self.files = [(0, PREFIX.stat().st_size, PREFIX)]
        for path in DIRECTORY.glob('*-*.bin'):
            start = int(path.name.split('-')[0])
            self.files.append((start, start + path.stat().st_size, path))

    def readable(self): return True
    def seekable(self): return True
    def tell(self): return self.pos
    def seek(self, n, whence=0):
        self.pos = n if whence == 0 else self.pos + n if whence == 1 else SIZE + n
        return self.pos

    def read(self, n=-1):
        end = min(SIZE, self.pos + n) if n >= 0 else SIZE
        output = []
        while self.pos < end:
            choices = [(start, stop, path) for start, stop, path in self.files if start <= self.pos < stop]
            if not choices:
                raise MissingRange((self.pos, end))
            start, stop, path = max(choices, key=lambda v: v[1])
            count = min(end, stop) - self.pos
            with path.open('rb') as handle:
                handle.seek(self.pos - start)
                part = handle.read(count)
            if len(part) != count:
                raise ValueError('Cached range truncated during read')
            self.pos += count
            output.append(part)
        return b''.join(output)

    def readinto(self, b):
        value = self.read(len(b)); b[:len(value)] = value; return len(value)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--extract', action='store_true')
    args = parser.parse_args()
    try:
        with h5py.File(CachedFile(), 'r') as h:
            print('obs keys', list(h['obs']), flush=True)
            for name in ['obs/gene', 'obs/gem_group', 'obs/__categories/gene', 'obs/__categories/gem_group', 'var/gene_name']:
                if name not in h: continue
                x = h[name]
                if isinstance(x, h5py.Dataset):
                    print(name, x.shape, str(x.dtype), x.id.get_offset(), x.chunks, x.compression, flush=True)
                else:
                    print(name, 'group', list(x), flush=True)
            if args.extract:
                from anndata.io import read_elem
                labels = read_elem(h['obs/gene'])
                batches = read_elem(h['obs/gem_group'])
                if 'obs/__categories' in h:
                    for key, value in [('gene', labels), ('gem_group', batches)]:
                        if key in h['obs/__categories']:
                            decoded = read_elem(h['obs/__categories'][key])[value]
                            if key == 'gene': labels = decoded
                            else: batches = decoded
                genes = read_elem(h['var/gene_name'])
                if 'var/__categories/gene_name' in h:
                    genes = read_elem(h['var/__categories/gene_name'])[genes]
                ds = h['X']
                assert ds.chunks is None and ds.compression is None and ds.dtype == np.float32
                np.savez_compressed(DIRECTORY/'metadata.npz', labels=np.asarray(labels,dtype=str),
                    batches=np.asarray(batches,dtype=str), genes=np.asarray(genes,dtype=str),
                    x_shape=np.asarray(ds.shape), x_offset=np.asarray(ds.id.get_offset()))
                print('SAVED metadata.npz', len(labels), len(genes), flush=True)
    except MissingRange as error:
        a, b = error.args[0]
        start = a // 1048576 * 1048576
        end = min(SIZE, max(start + 1048576, b))
        value = dict(start=start, end=end, requested_start=a, requested_end=b)
        (DIRECTORY/'missing.json').write_text(json.dumps(value))
        print('MISSING', json.dumps(value), flush=True)
        raise SystemExit(2)


if __name__ == '__main__': main()
