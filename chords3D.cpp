#include <complex>
#include <cstddef>
#include <cstdio>
#include <malloc/_malloc.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/_types/_null.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <omp.h>

#include "common.h"

struct chord_pixel {
    chord_pixel* previous_pixel;
    unsigned long pixel_index;
};

void set_chord(unsigned long* chords, chord_pixel last_pixel, unsigned int chord_length);
void chords_along_line(volume vol, unsigned long* chords, point line_start, point line_end);
unsigned long* chords_in_any_direction(volume vol, line l);
void merge_chords(unsigned long* chords, unsigned long* to_merge, unsigned int num_elements);
int save_to_file(unsigned long* chords, unsigned int num_elements);

int main() {
    // Dimensions of the 3D array
    int width = 500;  // replace with actual width
    int height = 500; // replace with actual height
    int depth = 500;  // replace with actual depth
    unsigned short number_of_angles = 1000;

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

    printf("Beginning processing\n");
    volume vol = {(unsigned short) width,(unsigned short) height,(unsigned short) depth, data };
    point center = {(unsigned short) (width / 2),(unsigned short) (height / 2),(unsigned short) (depth / 2)};

    unsigned long* all_chords = (unsigned long*) malloc((unsigned int) vol.width * vol.height * vol.depth * sizeof(unsigned long));

    // 1. Calculate fibonacci Sphere Points
    point fibonacci_points[number_of_angles];
    unsigned short fibonacci_radius = min(width, height, depth) / 2;
    fibonacci_sphere(fibonacci_points, number_of_angles, center, fibonacci_radius);

    // For Every Sphere Point
    # pragma omp parallel for
    for (int i = 0; i < number_of_angles; i++) {
        printf("Working on angle no.: %d\n", i);


        // 2. Create a line, by specifying generated point and center of sphere/volume
        line l = {center, fibonacci_points[i]};

        // printf("Processing fibonacci point: %d, %d, %d, %d, %d, %d\n", l.start.x, l.start.y, l.start.z, l.end.x, l.end.y, l.end.z);

        unsigned long* result = chords_in_any_direction(vol, l);


        merge_chords(all_chords, result, num_elements);
        free(result);

        if (i % 100 == 0) {
            printf("Saving to file\n");
            save_to_file(all_chords, num_elements);
        }
    }

    // save to file final result
    save_to_file(all_chords, num_elements);

    // Free the memory
    free(all_chords);

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

/*
 Generate chords along a specified lines direction
 Generate all parallel lines for this and then calculate chord length for each

 CURRENTLY ONLY FILLS ENTIRE VOLUME, IF LINE GOES THROUGH CENTER OF VOLUME
 AND ONLY WORKS WITH EQUAL X,Y,Z DIMENSIONS OF VOLUME
*/
unsigned long* chords_in_any_direction(volume vol, line l) {
    unsigned long* chords = (unsigned long*) malloc((unsigned int) vol.width * vol.height * vol.depth * sizeof(unsigned short));

    // Check if memory allocation was successful
    if (chords == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    // printf("Allocated memory\n");

    parametric_line pl = to_parametric_line(l);

    // printf("Converted to parametric line\n");

    unsigned int grid_size = vol.width;
    parametric_line* parallel_lines = get_parallel_lines(pl, grid_size);

    unsigned long num_paralell_lines = (grid_size * 2) * (grid_size * 2);

    // printf("Got parallel lines\n");

    for (int k = 0; k <= num_paralell_lines; k++) {

        double t_min, t_max;
        get_line_bounds(parallel_lines[k], vol, &t_min, &t_max);
        if (t_min == -1) {
            continue;
        }

        // Calculate the start and stop points of the chord
        point start = {(unsigned short) (parallel_lines[k].x0 + t_min * parallel_lines[k].direction.a),
                       (unsigned short) (parallel_lines[k].y0 + t_min * parallel_lines[k].direction.b),
                       (unsigned short) (parallel_lines[k].z0 + t_min * parallel_lines[k].direction.c)};
        point stop = {(unsigned short) (parallel_lines[k].x0 + t_max * parallel_lines[k].direction.a),
                      (unsigned short) (parallel_lines[k].y0 + t_max * parallel_lines[k].direction.b),
                      (unsigned short) (parallel_lines[k].z0 + t_max * parallel_lines[k].direction.c)};

        // DEBUG
        if (stop.x > 499 || stop.y > 499 || stop.z > 499 || start.x > 499 || start.y > 499 || start.z > 499) {
            //printf("Tried getting chords along line: %d, %d, %d, %d, %d, %d\n", start.x, start.y, start.z, stop.x, stop.y, stop.z);
            continue;
        }

        chords_along_line(vol, chords, start, stop);
    }

    free(parallel_lines);

    return chords;
}


// adjusted Bresenham's line algorithm to generate chords along a line and measure their length
void chords_along_line(volume vol, unsigned long* chords, point line_start, point line_end) {
    int x0 = line_start.x;
    int y0 = line_start.y;
    int z0 = line_start.z;
    int x1 = line_end.x;
    int y1 = line_end.y;
    int z1 = line_end.z;

    int dx = abs(x1-x0), sx = x0<x1 ? 1 : -1;
    int dy = abs(y1-y0), sy = y0<y1 ? 1 : -1;
    int dz = abs(z1-z0), sz = z0<z1 ? 1 : -1;
    int dm = max(dx,dy,dz), i = dm; /* maximum difference */
    x1 = y1 = z1 = dm/2; /* error offset */

    chord_pixel* last_chord_pixel = 0;
    unsigned int chord_length = 0;

    while (1) {  /* loop */
        unsigned int index = z0 * vol.width * vol.height + x0 * vol.height + y0;
        //printf("index: %d, %d, %d, %d, %d, %d, %d\n", index, z0, x0, y0, vol.height, vol.width, vol.depth);
        unsigned short value = vol.data[index];
        //printf("value: %d\n", value);

        if (value == 0) {
                //printf("Value 0 \n");
                chord_pixel* new_chord_pixel = (chord_pixel*) malloc(sizeof(chord_pixel));
                new_chord_pixel->previous_pixel = last_chord_pixel;
                new_chord_pixel->pixel_index = index;
                last_chord_pixel = new_chord_pixel;
        } else {
            if (last_chord_pixel != 0) {
                set_chord(chords, *last_chord_pixel, chord_length);

                // Free memory
                while (1) {
                    chord_pixel* to_free = last_chord_pixel;
                    if (last_chord_pixel->previous_pixel != 0) {
                        last_chord_pixel = last_chord_pixel->previous_pixel;
                        free(to_free);
                    } else {
                        break;
                        free(to_free);
                    }
                }
                chord_length = 0;
            }
        }
        chord_length++;

        if (i-- == 0) break;
        x1 -= dx; if (x1 < 0) { x1 += dm; x0 += sx; }
        y1 -= dy; if (y1 < 0) { y1 += dm; y0 += sy; }
        z1 -= dz; if (z1 < 0) { z1 += dm; z0 += sz; }
    }

    if (last_chord_pixel != 0) {
        set_chord(chords, *last_chord_pixel, chord_length);
    }

}

// use a stack to store chord pixels and set them efficiently
void set_chord(unsigned long* chords, chord_pixel last_pixel, unsigned int chord_length) {
    // Free up memory on the go, except for the last pixel, that will be freed by the caller
    while (last_pixel.previous_pixel != 0) {
        last_pixel = *last_pixel.previous_pixel;
        chords[last_pixel.pixel_index] = chord_length;
    }
}

void merge_chords(unsigned long* chords, unsigned long* to_merge, unsigned int num_elements) {
    for (int i = 0; i < num_elements; i++) {
        #pragma omp atomic
        chords[i] = chords[i] + to_merge[i];
    }
}

int save_to_file(unsigned long* chords, unsigned int num_elements) {
    // save to file
    FILE *fp = fopen("chords.raw", "wb");
    if (fp == NULL) {
        perror("Error opening file");
        return 1;
    }
    fwrite(chords, sizeof(unsigned long), num_elements, fp);
    fclose(fp);

    return 0;
}
