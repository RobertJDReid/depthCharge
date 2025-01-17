
# DepthCharge v0.18

![sub](sub2.png)

**DepthCharge** is a command-line tool written in C for calculating the mean sequence depth of genomic data within user-defined bins. It processes the output from `samtools depth` and supports the use of FASTA index files for chromosome length reference.

## Features

- Calculates mean sequence depth for a specified bin size (default: 1000).
- Handles chromosome transitions seamlessly.
- Supports optional FASTA index (`.fai`) files for accurate chromosome length handling.
- Processes large genomic datasets efficiently.

## Requirements

- C compiler (e.g., `gcc`)
- Genomic depth data generated by `samtools depth`
- Optional: A FASTA index (`.fai`) file for more accurate chromosome length determination.

## Installation

Compile the program using a C compiler:
```bash
gcc -o depthCharge depthCharge18.c -lm
```

link into `usr/local/bin` or a folder in your `$PATH` with `ln -s _path-to-depthCharge_`

## Usage

```bash
depthCharge [OPTIONS]
```

### Options

- `--bin-size=<integer>`  
  Set the size of the bins for calculating mean depth (default: 1000).

- `-i <file>`  
  Specify the input file. If not provided, the program reads from `stdin`.

- `--fai <file>` or `--fai=<file>`  
  Provide the path to a `.fai` file to reference chromosome lengths.

### Example Commands

1. **Basic Usage**  
   Calculate mean depth with the default bin size (1000):
   ```bash
   samtools depth input.bam | depthCharge
   ```

1. **Target to bed file**
   ```bash
   samtools depth input.bam | depthCharge > binnedDepth.bed
   ```

2. **Custom Bin Size**  
   Calculate mean depth with a bin size of 500:
   ```bash
   samtools depth input.bam | depthCharge --bin-size=500
   ```

3. **Using a FASTA Index File**  
   Calculate mean depth with a FASTA index file:
   ```bash
   samtools depth input.bam | depthCharge --fai genome.fai
   ```

4. **Using an Input File**  
   Process an existing depth file:
   ```bash
   depthCharge -i depth.txt --bin-size=1000 --fai genome.fai
   ```

## Input Format

The input should be in the format produced by `samtools depth`:
```
<chromosome>  <position>  <depth>
```

## Output Format

The output is tab-delimited and includes:
```
<chromosome>  <bin_start>  <bin_end>  <mean_depth>
```

## Notes

- Ensure that the input file and the `.fai` file (if used) are compatible and contain the required data.
- Chromosome transitions are handled, ensuring bins at the ends of chromosomes are correctly processed.

## License

This project is licensed under the MIT License.

## Author

**Robert J.D. Reid**  
Feel free to contact me for questions or feature requests
