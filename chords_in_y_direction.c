#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

void set_chord(unsigned short* mask, unsigned int chord_start, unsigned int chord_stop) {
    for (unsigned int i = chord_start; i < chord_stop; i++) {
        mask[i] = chord_stop - chord_start + 1;
    }
}

unsigned short* chords_in_y_direction(unsigned char* mask, int width, int height, int depth)
{

    unsigned short* chords = (unsigned short*) malloc((unsigned int) height * depth * width * sizeof(unsigned short));

    // Check if memory allocation was successful
        if (chords == NULL) {
            fprintf(stderr, "Memory allocation failed\n");
            return NULL;
        }

    printf("Beginning processing\n");

    bool in_chord = false;
    unsigned int chord_start = 0;

    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < width; j++) {
            for (int k= 0; k < height; k++) {
                unsigned int index = i * width * height + j * height + k;
                if (mask[index] == 0) {
                    if (!in_chord) {
                        in_chord = true;
                        chord_start = index;
                    }
                    continue;
                } else {
                    if (in_chord) {
                        in_chord = false;
                        set_chord(chords, chord_start, index - 1);
                    }
                }
            }
            if (in_chord) {
                in_chord = false;
                set_chord(chords, chord_start, (unsigned int) (i * width * height + (j + 1) * height) - 1);
            }
        }
        if (in_chord) {
            in_chord = false;
            set_chord(chords, chord_start, (unsigned int) ((i + 1) * width * height) - 1);
        }
    }

    if (in_chord) {
        in_chord = false;
        set_chord(chords, chord_start, (unsigned int) depth * width * height - 1);
    }

    printf("Chords layed through volume! \n");

    return chords;
}

int main() {
    // Dimensions of the 3D array
    int width = 500;  // replace with actual width
    int height = 500; // replace with actual height
    int depth = 500;  // replace with actual depth

    // Calculate the total number of elements
    unsigned int num_elements = width * height * depth;
    unsigned int file_size = num_elements;

    // Open the binary file
    int fd = open("../lung424_small_mask.raw", O_RDONLY);
    if (fd == -1) {
        perror("Error opening file");
        return 1;
    }

   	struct stat file_properties;
	if(fstat(fd, &file_properties) < 0)
	{
		perror("Error getting file properties");
		return 1;
	}

	unsigned int N_bytes = file_properties.st_size;
	printf("File size: %u\n", (unsigned int) N_bytes);
	printf("Expected size: %u\n", (unsigned int) file_size);

	if ((unsigned int) N_bytes != (unsigned int) file_size)
    {
        printf("File size does not match the expected size\n");
        return 1;
    }

    // Memory map the file
    unsigned char *data = (unsigned char *)mmap(NULL, (unsigned int) file_size, PROT_READ, MAP_SHARED, fd, 0);
    if (data == MAP_FAILED) {
        perror("Error mmapping the file");
        close(fd);
        return 1;
    }

    unsigned short* result = chords_in_y_direction(data, width, height, depth);

    // save to file
    FILE *fp = fopen("chords.raw", "wb");
    if (fp == NULL) {
        perror("Error opening file");
        return 1;
    }
    for (int i = 0; i < depth; i++) {
        for (int j = 0; j < width; j++) {
            for (int k = 0; k < height; k++) {
                unsigned int index = i * width * height + j * height + k;
                fwrite(&result[index], sizeof(unsigned short), 1, fp);
            }
        }
    }
    fclose(fp);

    // Free the memory
    free(result);

    // Unmap the file
    if (munmap(data, file_size) == -1) {
        perror("Error unmapping the file");
        close(fd);
        return 1;
    }

    // Close the file
    close(fd);

    return 0;
}
