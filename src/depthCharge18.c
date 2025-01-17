#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

// **************** depthCharge version 0.18
// Takes a file of sequence depth by position output from samtoosl depth and returns
// mean sequence depth for a user-defined bin size.
// Default bin == 1000, but can be set arbitrarily with --bin-size=integer
//
// fixes chromosome change so last bin with data is properly handled (not zeroed out before processing)
// EOF fixed
// loading of chromosome sizes from index fixed

#define DEFAULT_BIN_SIZE 1000
#define MAX_CHROMOSOMES 1000

// Structure to store chromosome lengths
typedef struct {
    char name[256];
    int length;
} ChromosomeInfo;

// Function to calculate mean coverage for a bin and print it
void process_bin(const char *chromosome, int start, int end, long long total_coverage, int bin_size) {
    if (chromosome[0] != '\0') { // Ensure chromosome is valid
        double mean_coverage = (double) total_coverage / bin_size;
        printf("%s\t%d\t%d\t%.2f\n", chromosome, start, end, mean_coverage);
    }
}

// Function to read the FASTA index file and load chromosome lengths
int load_fasta_index(const char *fai_path, ChromosomeInfo *chromosomes) {
    FILE *fai_file = fopen(fai_path, "r");
    if (fai_file == NULL) {
        fprintf(stderr, "Error: Unable to open FASTA index file %s.\n", fai_path);
        return -1;
    }

    int count = 0;
    char name[256];
    int length;
    // Variables for unused columns
    long offset;
    int line_bases, line_width;

    while (fscanf(fai_file, "%255s\t%d\t%ld\t%d\t%d", name, &length, &offset, &line_bases, &line_width) == 5) {
        // Copy chromosome name and length into the data structure
        strncpy(chromosomes[count].name, name, sizeof(chromosomes[count].name) - 1);
        chromosomes[count].name[sizeof(chromosomes[count].name) - 1] = '\0';
        chromosomes[count].length = length;

        // Debug: Print each chromosome name and length as it's loaded
        //fprintf(stderr, "DEBUG: Loaded chromosome %s with length %d\n", chromosomes[count].name, chromosomes[count].length);

        count++;
        if (count >= MAX_CHROMOSOMES) {
            fprintf(stderr, "Error: Exceeded maximum number of chromosomes (%d).\n", MAX_CHROMOSOMES);
            fclose(fai_file);
            return -1;
        }
    }

    fclose(fai_file);
    //fprintf(stderr, "FAI file loaded successfully with %d chromosomes.\n", count);
    return count;
}


// Function to get the length of a chromosome from the loaded index
int get_chromosome_length(const ChromosomeInfo *chromosomes, int chromosome_count, const char *chromosome) {
    for (int i = 0; i < chromosome_count; i++) {
        if (strcmp(chromosomes[i].name, chromosome) == 0) {
            // DEBUG
            //fprintf(stderr, "Current chromosome length from fai %d\n",chromosomes[i].length);
            return chromosomes[i].length;
        }
    }
    return -1; // Chromosome not found
}

int main(int argc, char *argv[]) {
    int bin_size = DEFAULT_BIN_SIZE;
    FILE *input_file = stdin;
    char fai_path[256] = "";
    bool use_fai = false;

    // Parse arguments for --bin-size, -i, and --fai

for (int i = 1; i < argc; i++) {
    if (strncmp(argv[i], "--bin-size=", 11) == 0) {
        bin_size = atoi(argv[i] + 11);
        if (bin_size <= 0) {
            fprintf(stderr, "Error: Invalid bin size. Must be a positive integer.\n");
            return EXIT_FAILURE;
        }
    } else if (strcmp(argv[i], "-i") == 0) {
        if (i + 1 < argc) {
            input_file = fopen(argv[i + 1], "r");
            if (input_file == NULL) {
                fprintf(stderr, "Error: Unable to open file %s.\n", argv[i + 1]);
                return EXIT_FAILURE;
            }
            i++;
        } else {
            fprintf(stderr, "Error: -i option requires a file name.\n");
            return EXIT_FAILURE;
        }
    } else if (strcmp(argv[i], "--fai") == 0) {
        if (i + 1 < argc) {
            strncpy(fai_path, argv[i + 1], sizeof(fai_path) - 1);
            fai_path[sizeof(fai_path) - 1] = '\0';
            use_fai = true;
            i++; // Skip the next argument as it's the path
        } else {
            fprintf(stderr, "Error: --fai option requires a file path.\n");
            return EXIT_FAILURE;
        }
    } else if (strncmp(argv[i], "--fai=", 6) == 0) {
        strncpy(fai_path, argv[i] + 6, sizeof(fai_path) - 1);
        fai_path[sizeof(fai_path) - 1] = '\0';
        use_fai = true;
    }
}


    ChromosomeInfo chromosomes[MAX_CHROMOSOMES];
    int chromosome_count = 0;
    if (use_fai) {
        chromosome_count = load_fasta_index(fai_path, chromosomes);
        if (chromosome_count < 0) {
            return EXIT_FAILURE;
        }
    }

    char current_chromosome[256] = "";
    int current_bin_start = 1; // 1-based bins
    int current_bin_end = bin_size;
    long long total_coverage = 0;

    char chromosome[256];
    int position;
    int coverage;
    int last_position = 0;

    while (fscanf(input_file, "%s\t%d\t%d", chromosome, &position, &coverage) == 3) {
        // DEBUG: Print the current line being processed
        // fprintf(stderr, "DEBUG: Processing line: %s\t%d\t%d\n", chromosome, position, coverage);
    
        // Check if first run where current_chromosome is null
        if (current_chromosome[0] == '\0') {
            strncpy(current_chromosome, chromosome, sizeof(current_chromosome) - 1);
            current_chromosome[sizeof(current_chromosome) - 1] = '\0';
        }

        // Handle bin transitions and chromosome changes
        while (strcmp(current_chromosome, chromosome) != 0 || position > current_bin_end) {
            if (strcmp(current_chromosome, chromosome) != 0 && current_chromosome[0] != '\0') {
                int chromosome_length = use_fai ? get_chromosome_length(chromosomes, chromosome_count, current_chromosome) : last_position;
                if (chromosome_length < 0) {
                    chromosome_length = last_position; // Default to last position if length is unknown
                    //fprintf(stderr, "Warning: Using last position (%d) as fallback for chromosome %s.\n", last_position, current_chromosome);
                }
                // DEBUG: Print the current state of the loop
                //fprintf(stderr, "DEBUG: Current chromosome: %s, Position: %d, Bin: [%d-%d], Chromosome Length: %d, Total: %lld\n",
                //    current_chromosome, position, current_bin_start, current_bin_end, chromosome_length, total_coverage);

                while (current_bin_start <= chromosome_length) {
                    int partial_bin_end = (current_bin_end > chromosome_length) ? chromosome_length : current_bin_end;
                    int partial_bin_size = (partial_bin_end - current_bin_start + 1);
                    process_bin(current_chromosome, current_bin_start, partial_bin_end, total_coverage, partial_bin_size); // problem here
                    current_bin_start = current_bin_end + 1;
                    current_bin_end = current_bin_start + bin_size - 1;
                    total_coverage = 0;
                }

                strncpy(current_chromosome, chromosome, sizeof(current_chromosome) - 1);
                current_chromosome[sizeof(current_chromosome) - 1] = '\0';
                current_bin_start = 1;
                current_bin_end = bin_size;
                total_coverage = 0;
                break;
            }

            process_bin(current_chromosome, current_bin_start, current_bin_end, total_coverage, bin_size);
            current_bin_start = current_bin_end + 1;
            current_bin_end = current_bin_start + bin_size - 1;
            total_coverage = 0;
        }

        total_coverage += coverage;
        last_position = position;
    }

    if (current_chromosome[0] != '\0') {
        int chromosome_length = use_fai ? get_chromosome_length(chromosomes, chromosome_count, current_chromosome) : last_position;
        if (chromosome_length < 0) chromosome_length = last_position;
        while (current_bin_start <= chromosome_length) {
            int partial_bin_end = (current_bin_end > chromosome_length) ? chromosome_length : current_bin_end;
            int partial_bin_size = (partial_bin_end - current_bin_start + 1);
            process_bin(current_chromosome, current_bin_start, partial_bin_end, total_coverage, partial_bin_size);
            current_bin_start = current_bin_end + 1;
            current_bin_end = current_bin_start + bin_size - 1;
            total_coverage = 0;
        }
    }

    if (input_file != stdin) {
        fclose(input_file);
    }

    return EXIT_SUCCESS;
}
