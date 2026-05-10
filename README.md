# Detection of local motifs in DNA sequences

This project implements a polynomial-time heuristic algorithm for finding common motifs in nucleotide sequences using graph-based methods.

The algorithm:

- reads nucleotide sequences with quality/confidence values,
- removes positions with confidence below a user-defined threshold,
- generates all k-mers (substrings of length `4–9`),
- builds a graph where:
  - each k-mer occurrence is a separate vertex,
  - edges connect identical k-mers from different sequences if their positions are sufficiently close,
- heuristically searches for a clique-like structure containing exactly one vertex from each sequence.

The output contains:

- sequence identifier,
- original position in the input sequence,
- detected nucleotide substring.

The program was tested with different confidence thresholds and k-mer lengths to analyze their impact on graph structure and motif detection results.

The heuristic is fully graph-based and does not use biological sequence information during the clique search stage.

## Usage

```bash
c++ 'AKwB 2.cpp' -o ./a.out
./a.out
```
