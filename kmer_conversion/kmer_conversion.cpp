#include <iostream>
#include <string>
#include <chrono>
// A function that gets a kmer, in the form of a string, and transforms iot into a 2bit number

typedef struct {
    unsigned int value : 2; // 2-bit wide field
} EncodedBase;


EncodedBase encode_base(char base) {
    EncodedBase encoded_base;
    // A in ASCII is:           01000  00  1
    // C in ASCII is:           01000  01  1
    // G in ASCII is:           01000  11  1
    // T in ASCII is:           01010  10  0
    encoded_base.value = (base >> 1) % 4;
    return encoded_base;
}


void encode_read(const char* read, size_t length, EncodedBase* encoded_read) {
    for (int i = 0; i < length; i++) encoded_read[i].value = encode_base(read[i]).value;
}

uint64_t* read_file(const char *filePath) {
    FILE *file = fopen(filePath, "r"); // Open the file for reading
    if (file == NULL) {
        perror("Error opening file");
        return 0;
    }

    char *line = NULL;
    size_t len = 0;
    ssize_t read;
    uint64_t result = 0;
    size_t kmer_size = 32;
    uint64_t base_kmer_addr = 0;
    uint64_t* addresses_to_search = new uint64_t[kmer_size * 2];

    uint64_t precomputed_offsets[64];
    for (size_t j = 0; j < 64; j++) {
        precomputed_offsets[j] = 2 * j;
    }


    EncodedBase*  encoded_read = new EncodedBase[1000];
    size_t line_number = 0;
    while ((read = getline(&line, &len, file)) != -1) {
        if (line_number++ % 4 != 1) {
            continue; // Skip all lines that are not reads
        }
        encode_read(line, len, encoded_read);
        size_t num_kmers = len - kmer_size + 1;
        // Assuming you have a precomputed table for offsets
    // uint64_t precomputed_offsets[64];

    for (size_t i = 0; i < num_kmers; i++) {
        uint64_t kmer_value = 0;
        for (size_t j = 0; j < kmer_size; j += 4) {
            kmer_value = (kmer_value << 8) | 
                        (encoded_read[i + j].value << 6) |
                        (encoded_read[i + j + 1].value << 4) |
                        (encoded_read[i + j + 2].value << 2) |
                        encoded_read[i + j + 3].value;
        }

        for (size_t j = 0; j < 64; j++) {
            addresses_to_search[j] += base_kmer_addr + precomputed_offsets[j] + ((kmer_value >> j) & 1);
        }
    }

    }
    delete[] encoded_read;

    free(line); // Free the allocated buffer
    fclose(file); // Close the file
    return addresses_to_search;
}


int main() {
    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    // Call the function you want to time
    read_file("../data/reads/illumina-sars-cov-2.fq");

    // Stop the timer
    auto stop = std::chrono::high_resolution_clock::now();

    // Calculate the duration
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Time taken by function: " << duration.count() << " milliseconds" << std::endl;

    return 0;
}