# Dashing2 Python Bindings

Python bindings for the Dashing2 sequence sketching library.

## Installation

To install the Dashing2 Python bindings, run the following command:

```bash
pip install dashing2
```

## Usage

To use the Dashing2 Python bindings, import the library and use the provided classes and functions.

```python
from dashing2 import RMinHash, HLL, sketch_fastx
# Using MinHash Sketching
mh = RMinHash(sketch_size=1024)
mh.add_sequence("ACGTACGT", kmer_size=3)

# Using HyperLogLog Sketching
hll = HLL(precision=14)
hll.add_sequence("ACGTACGT", kmer_size=3)
estimate = hll.get_estimate()

# Sketching from FASTA/FASTQ files
result = sketch_fastx("path/to/sequences.fa")
```

## Features

- MinHash sketching with RMinHash
- HyperLogLog cardinality estimation
- FASTA/FASTQ file processing
- Efficient k-mer based sequence comparison
- Sketch merging and comparison operations

## Requirements

- Python >= 3.6
- numpy >= 1.15.0