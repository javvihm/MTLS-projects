# Dataset

This dataset contains 1,200 protein entries, each provided in three complementary file formats within dedicated subfolders. All entries share the same base filename across subfolders.

The dataset structure is as follows:

```bash
dataset/
├── dssp/
│   ├── 1a00.dssp
│   ├── 1a01.dssp
│   ├── ...
├── fasta/
│   ├── 1a00.fasta
│   ├── 1a01.fasta
│   ├── ...
└── pssm/
    ├── 1a00.pssm
    ├── 1a01.pssm
    ├── ...
```

## Subfolders and File Formats

### dssp

Contains files in DSSP format, which hold secondary structure assignments computed by the DSSP algorithm. Each file lists secondary structure elements (e.g., H for helix, E for strand, and C for coil).

### fasta

Contains files in FASTA format, the standard format for representing biological sequences. Each file starts with a header (beginning with >) followed by the protein sequence.

### pssm

Contains files in PSSM (Position-Specific Scoring Matrix) format. These CSV files represent evolutionary profiles for proteins with a matrix of weighted observed percentages for each residue position (rounded down).
