#!/usr/bin/env python3
"""
dna_analysis.py - A simple bioinformatics tool for DNA sequence analysis.

This script provides functions to calculate GC content, reverse complements,
motif finding, melting temperature, and other properties for DNA sequences from
FASTA files with multiple sequences. Includes error handling for invalid bases.

Author: [Swarnendu Das]
Date: April 04, 2025
"""

def validate_dna_sequence(dna_sequence):
    """
    Validate that a DNA sequence contains only A, T, G, C.

    Args:
        dna_sequence (str): Input DNA sequence

    Raises:
        ValueError: If the sequence contains invalid bases
    """
    valid_bases = set('ATGC')
    sequence_upper = dna_sequence.upper()
    invalid_bases = set(sequence_upper) - valid_bases
    
    if invalid_bases:
        raise ValueError(f"Invalid bases found: {invalid_bases}. Only A, T, G, C are allowed.")

def read_fasta_file(file_path):
    """
    Read multiple DNA sequences from a FASTA file.

    Args:
        file_path (str): Path to the FASTA file

    Returns:
        dict: Dictionary with header (without '>') as key and sequence as value

    Raises:
        FileNotFoundError: If the file does not exist
        ValueError: If the file is empty or has no sequence data
    """
    try:
        with open(file_path, 'r') as file:
            sequences = {}
            current_header = None
            current_sequence = ''
            
            for line in file:
                line = line.strip()
                if line.startswith('>'):
                    if current_header and current_sequence:
                        sequences[current_header] = current_sequence
                    current_header = line[1:]  # Remove '>' from header
                    current_sequence = ''
                elif line:
                    current_sequence += line
            
            if current_header and current_sequence:
                sequences[current_header] = current_sequence
            
            if not sequences:
                raise ValueError("No sequence data found in FASTA file.")
            return sequences
    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file not found: {file_path}")

def calculate_gc_content(dna_sequence):
    """
    Calculate the GC content percentage of a DNA sequence.

    Args:
        dna_sequence (str): Input DNA sequence (e.g., "ATGC")

    Returns:
        float: GC percentage, or 0 if sequence is empty

    Raises:
        ValueError: If the sequence contains invalid bases
    """
    if not dna_sequence:
        return 0
    
    validate_dna_sequence(dna_sequence)
    dna_sequence = dna_sequence.upper()
    gc_count = dna_sequence.count('G') + dna_sequence.count('C')
    total_length = len(dna_sequence)
    return (gc_count / total_length) * 100

def calculate_at_content(dna_sequence):
    """
    Calculate the AT content percentage of a DNA sequence.

    Args:
        dna_sequence (str): Input DNA sequence (e.g., "ATGC")

    Returns:
        float: AT percentage, or 0 if sequence is empty

    Raises:
        ValueError: If the sequence contains invalid bases
    """
    if not dna_sequence:
        return 0
    
    validate_dna_sequence(dna_sequence)
    dna_sequence = dna_sequence.upper()
    at_count = dna_sequence.count('A') + dna_sequence.count('T')
    total_length = len(dna_sequence)
    return (at_count / total_length) * 100

def calculate_melting_temp(dna_sequence):
    """
    Calculate the melting temperature (Tm) of a DNA sequence using the Wallace rule.
    Tm = 2°C * (A+T) + 4°C * (G+C). Suitable for short sequences (< 50 bp).

    Args:
        dna_sequence (str): Input DNA sequence (e.g., "ATGC")

    Returns:
        float: Melting temperature in °C, or 0 if sequence is empty

    Raises:
        ValueError: If the sequence contains invalid bases
    """
    if not dna_sequence:
        return 0
    
    validate_dna_sequence(dna_sequence)
    dna_sequence = dna_sequence.upper()
    at_count = dna_sequence.count('A') + dna_sequence.count('T')
    gc_count = dna_sequence.count('G') + dna_sequence.count('C')
    return (2 * at_count) + (4 * gc_count)

def calculate_molecular_weight(dna_sequence):
    """
    Calculate the approximate molecular weight of a single-stranded DNA sequence.
    Uses average weights: A=313.2, T=304.2, G=329.2, C=289.2 g/mol.

    Args:
        dna_sequence (str): Input DNA sequence (e.g., "ATGC")

    Returns:
        float: Molecular weight in g/mol, or 0 if sequence is empty

    Raises:
        ValueError: If the sequence contains invalid bases
    """
    if not dna_sequence:
        return 0
    
    validate_dna_sequence(dna_sequence)
    dna_sequence = dna_sequence.upper()
    weights = {'A': 313.2, 'T': 304.2, 'G': 329.2, 'C': 289.2}
    return sum(weights[base] for base in dna_sequence)

def reverse_complement(dna_sequence):
    """
    Generate the reverse complement of a DNA sequence.

    Args:
        dna_sequence (str): Input DNA sequence (e.g., "ATGC")

    Returns:
        str: Reverse complement sequence

    Raises:
        ValueError: If the sequence contains invalid bases
    """
    validate_dna_sequence(dna_sequence)
    complements = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    dna_sequence = dna_sequence.upper()
    return ''.join(complements[base] for base in reversed(dna_sequence))

def find_motif(dna_sequence, motif):
    """
    Find all positions of a specific motif in the DNA sequence.

    Args:
        dna_sequence (str): Input DNA sequence
        motif (str): Short sequence to search for

    Returns:
        list: List of 1-based positions where motif occurs

    Raises:
        ValueError: If the sequence or motif contains invalid bases
    """
    validate_dna_sequence(dna_sequence)
    validate_dna_sequence(motif)
    
    positions = []
    motif_length = len(motif)
    dna_sequence = dna_sequence.upper()
    motif = motif.upper()
    
    for i in range(len(dna_sequence) - motif_length + 1):
        if dna_sequence[i:i + motif_length] == motif:
            positions.append(i + 1)  # 1-based indexing
    return positions

def analyze_sequence(dna_sequence, motif, header=None):
    """
    Perform full analysis on a DNA sequence and print results.

    Args:
        dna_sequence (str): DNA sequence to analyze
        motif (str): Motif to search for
        header (str, optional): Sequence identifier for display
    """
    gc_percentage = calculate_gc_content(dna_sequence)
    at_percentage = calculate_at_content(dna_sequence)
    melting_temp = calculate_melting_temp(dna_sequence)
    mol_weight = calculate_molecular_weight(dna_sequence)
    rev_comp = reverse_complement(dna_sequence)
    motif_positions = find_motif(dna_sequence, motif)
    base_counts = {
        'A': dna_sequence.count('A'),
        'T': dna_sequence.count('T'),
        'G': dna_sequence.count('G'),
        'C': dna_sequence.count('C')
    }
    
    prefix = f"Sequence '{header}': " if header else "Default Sequence: "
    print(f"{prefix}")
    print(f"  DNA Sequence: {dna_sequence}")
    print(f"  GC Content: {gc_percentage:.2f}%")
    print(f"  AT Content: {at_percentage:.2f}%")
    print(f"  Melting Temperature: {melting_temp:.1f}°C (Wallace rule)")
    print(f"  Molecular Weight: {mol_weight:.1f} g/mol")
    print(f"  Reverse Complement: {rev_comp}")
    print(f"  Positions of motif '{motif}': {motif_positions}")
    print(f"  Sequence Length: {len(dna_sequence)}")
    print(f"  Base Counts: {base_counts}\n")

def main():
    """Main function to demonstrate DNA sequence analysis with multi-sequence FASTA support."""
    # Default sample DNA sequence (used if no file is provided)
    sample_dna = "ATGGCCATAGCTAGCTAGCTAGCGCGCGCTA"
    motif = "GCTA"
    fasta_file = "sample_dna.fasta"  # Example FASTA file name
    
    try:
        # Attempt to read from FASTA file if it exists, otherwise use sample_dna
        try:
            sequences = read_fasta_file(fasta_file)
            print(f"Loaded {len(sequences)} sequence(s) from {fasta_file}")
            for header, sequence in sequences.items():
                analyze_sequence(sequence, motif, header)
        except FileNotFoundError:
            print(f"{fasta_file} not found. Using default sample sequence.")
            analyze_sequence(sample_dna, motif)
        
        # Test with invalid sequence
        invalid_dna = "ATGCX"
        print("Testing invalid sequence:")
        calculate_gc_content(invalid_dna)
        
    except (ValueError, FileNotFoundError) as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
