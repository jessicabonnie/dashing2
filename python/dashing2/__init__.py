from .core import RMinHash, HLL, sketch_fastx, compute_distances

__version__ = '0.0.1'

__all__ = ['RMinHash', 'HLL', 'sketch_fastx', 'compute_distances']

def sketch_file_and_compute_distances(filepath, sketch_size=1024, kmer_size=31):
    """
    Sketch sequences from a file and compute their distances
    
    Args:
        filepath: Path to FASTA/FASTQ file
        sketch_size: Size of each sketch
        kmer_size: Size of k-mers
        
    Returns:
        tuple: (names, distances) where names is a list of sequence names
        and distances is a numpy array of pairwise distances
    """
    result = sketch_fastx(filepath, sketch_size, kmer_size)
    distances = compute_distances(result["sketches"])
    return result["names"], distances

def sketch_sequences(sequences, sketch_size=1024, kmer_size=31):
    """
    Sketch a list of sequences
    
    Args:
        sequences: List of DNA sequences
        sketch_size: Size of sketch
        kmer_size: Size of k-mers
        
    Returns:
        numpy.ndarray: Array of sketches
    """
    sketches = np.zeros((len(sequences), sketch_size), dtype=np.uint64)
    
    for i, seq in enumerate(sequences):
        sketch = RMinHash(sketch_size)
        sketch.add_sequence(seq, kmer_size)
        sketches[i] = sketch.get_mins()
    
    return sketches
