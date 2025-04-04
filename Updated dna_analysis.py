#!/usr/bin/env python3
"""
dna_analysis.py - A simple bioinformatics tool for DNA sequence analysis.

This script provides functions to calculate GC content, generate reverse complements,
and find motifs in DNA sequences. Includes error handling for invalid bases.
Suitable for beginner-to-intermediate Python users interested in bioinformatics.

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

def main():
    """Main function to demonstrate DNA sequence analysis with error handling."""
    # Sample DNA sequence (valid)
    sample_dna = "ATGGCCATAGCTAGCTAGCTAGCGCGCGCTA"
    motif = "GCTA"
    
    try:
        # Perform analyses
        gc_percentage = calculate_gc_content(sample_dna)
        rev_comp = reverse_complement(sample_dna)
        motif_positions = find_motif(sample_dna, motif)
        base_counts = {
            'A': sample_dna.count('A'),
            'T': sample_dna.count('T'),
            'G': sample_dna.count('G'),
            'C': sample_dna.count('C')
        }
        
        # Display results
        print(f"DNA Sequence: {sample_dna}")
        print(f"GC Content: {gc_percentage:.2f}%")
        print(f"Reverse Complement: {rev_comp}")
        print(f"Positions of motif '{motif}': {motif_positions}")
        print(f"Sequence Length: {len(sample_dna)}")
        print(f"Base Counts: {base_counts}")
        
        # Test with invalid sequence
        invalid_dna = "ATGCX"
        print("\nTesting invalid sequence:")
        calculate_gc_content(invalid_dna)
        
    except ValueError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    main()
