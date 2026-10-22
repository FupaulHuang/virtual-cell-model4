# GWPS raw-count subset acquisition

Source: public Replogle K562 GWPS H5AD, Figshare file 35775507,
`https://ndownloader.figshare.com/files/35775507`.

The source object is 65,830,941,948 bytes. Its ETag is
`e85ededd1213eae9405a2f2422edef5d-252`. HDF5 inspection established an
uncompressed, contiguous float32 `X` with shape 1,989,578 × 8,248 and byte
offset 2,048. Metadata were read through a local file-like cache combining
the existing, read-only full-download prefix and explicit public byte ranges.
Gene and batch category codes were decoded before row selection.

The selected panel is the 575-target union of the 2026 challenge panel and
all 300 public 2025 H1 targets. With seed 20261020, selection takes up to
200 cells per observed target and up to 50 non-targeting cells per GEM batch.
The resulting 96,111 cells include 13,324 controls across 267 GEM batches.
There are 527 observed union targets (522 with at least 20 selected cells)
and 272 observed challenge targets (269 with at least 20 selected cells).
Targets absent from GWPS remain explicit zero-cell rows in the final bank.

`plan.json` describes 87,047 bounded ranges totaling 3,316,586,784 bytes.
Adjacent selected rows may share a range, with at most one intervening row
and a maximum range length of 2 MiB. Only the explicitly selected rows are
aggregated. Complete ranges already present in the original source prefix
were copied read-only into the range cache. The original full download was
neither interrupted nor overwritten by this workflow.

Figshare's anonymous S3 redirects expire after ten seconds. A preliminary
attempt to request a fresh Figshare redirect for every range encountered
HTTP 429 and was stopped. The production workflow requests one redirect
per batch of 256 S3 ranges, uses 24 parallel S3 requests, verifies each
range's exact expected length, and renews the redirect for the next batch.
Each request sends the source `If-Match` ETag. Expired requests are not
accepted as data; missing ranges remain pending under a fresh signature.

The network operation is always a direct approved `curl -L` invocation.
Python only creates configurations, reads cached metadata/counts, and
verifies artifacts. In this environment, adding shell output redirection
to the curl tool command caused sandbox DNS failures, so curl's own
`%output{>>...}` write-out facility records HTTP statuses instead.

Scripts:

- `scripts/gwps_range_metadata.py`: inspect cached HDF5 metadata and report
  missing byte ranges; `--extract` creates `metadata.npz`.
- `scripts/gwps_subset_bank.py plan`: select cells and generate the range plan.
- `scripts/gwps_signed_batch.py --start N --limit 256`: read freshly fetched
  `fresh_headers.txt` and emit `signed_batch.curl` without exposing the URL.
- `curl -L --parallel --parallel-max 24 --config data/raw/gwps_subset/signed_batch.curl`:
  retrieve one batch.
- `scripts/gwps_signed_batch.py --verify`: inspect exact file lengths and report
  the next starting index, retaining any missing ranges for retry.
- `scripts/gwps_subset_bank.py build`: aggregate the complete cache into the
  AtlasShift-compatible K562 bank and write its provenance receipt.

The production bank uses exact int64 count sums and arithmetic means of
per-cell CPM over all 8,248 measured source columns. Duplicate gene symbols,
if present, are summed. Batch-matched control probabilities use target
library-mass weights; batch-matched arithmetic CPM means use target cell-count
weights. Variances of the target and matched-control means are retained.

The legacy output name `K562_GWPS_CPM_full_statistics.npz` denotes a complete
target-axis statistics schema. It does not imply that every GWPS cell was
downloaded: the bank's `source` field and receipt explicitly identify the
selected-cell raw-count subset.
