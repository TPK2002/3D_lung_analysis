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
#include <chrono>
#include <iostream>
#include "common.h"

struct chord_pixel {
    chord_pixel* previous_pixel;
    unsigned long pixel_index;
};

void set_chord(unsigned long* chords, chord_pixel last_pixel, unsigned int chord_length);
void chords_along_line(volume vol, unsigned long* chords, point line_start, point line_end, bool chords_through, std::vector<unsigned short>& chord_lengths);
unsigned long* chords_in_any_direction(volume vol, line l, bool chords_through, std::vector<unsigned short>& chord_lengths);
void merge_chords(unsigned long* chords, unsigned long* to_merge, unsigned long num_elements, std::string mode);
unsigned short* compress_results(unsigned long* chords, unsigned long num_elements, unsigned short num_of_angles, std::string mode);

int main(int argc, char *argv[]) {
    // Get from CLI. nx, ny, nz, number_of_angles, mode
    std::string file_path;
    unsigned short nx, ny, nz, number_of_angles;
    std::string mode = "avg"; // default is average
    bool chords_through = 0; // should the processing happen on background: 0 or tissue: 1

    if (argc < 6) {
        printf("Chords 3D has to be called with arguments: <file_path> <nx> <ny> <nz> <num_angles> <mode> < --inverse>)");
        return 1;
    }

    // Parse CLI parameters
    // note: argv[0] is the programm name, thus begin with 1
    file_path = argv[1];
    nx = static_cast<unsigned short>(strtoul(argv[2],NULL,10));
    ny = static_cast<unsigned short>(strtoul(argv[3],NULL,10));
    nz = static_cast<unsigned short>(strtoul(argv[4],NULL,10));
    number_of_angles = static_cast<unsigned short>(strtoul(argv[5],NULL,10));

    if (argc >= 7) {
        mode = argv[6];
    }

    // if --inverse is specified as 7th parameter, invert processing => run processing on tissue
    if (argc >= 8) {
        if (strcmp("--inverse", argv[7]) == 0) {
            chords_through = 1;
        }
    }

    printf("Chords 3D called for: %s, %hux%hux%hu.\n", file_path.c_str(), nx, ny, nz);
    printf("Doing %hu angles per pixel with mode: %s \n",number_of_angles, mode.c_str());
    if (chords_through == 1) {
        printf("INVERSE MODE ACTIVE \n");
    }

    // Calculate the total number of elements
    unsigned long num_elements = nx * ny * nz;
    unsigned long file_size = num_elements;

    // Open the binary file
    int fd = open(file_path.c_str(), O_RDONLY);
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
	if ((unsigned int) N_bytes != (unsigned int) file_size)
    {
	    printf("File size: %u\n", (unsigned int) N_bytes);
	    printf("Expected size: %u\n", (unsigned int) file_size);
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
    volume vol = {(unsigned short) nx,(unsigned short) ny,(unsigned short) nz, data };
    point center = {(unsigned short) (nx / 2),(unsigned short) (ny / 2),(unsigned short) (nz / 2)};

    unsigned long* all_chords = (unsigned long*) malloc((unsigned int) vol.nx * vol.ny * vol.nz * sizeof(unsigned long));

    // If we use mode min, we need to set all_chords to SHRT_MAX
    if (mode == "min") {
        for (int i = 0; i < num_elements; i++) {
            all_chords[i] = ULONG_MAX;
        }
    }

    // 1. Calculate fibonacci Hemisphere Points
    point fibonacci_points[number_of_angles];
    unsigned short fibonacci_radius = min(nx, ny, nz) / 2;
    fibonacci_hemisphere(fibonacci_points, number_of_angles, center, fibonacci_radius);

    std::vector<std::vector<unsigned short>> loose_chord_lengths(number_of_angles);

    // For Every Hemisphere Point
    # pragma omp parallel for
    for (int i = 0; i < number_of_angles; i++) {
        printf("Working on angle no.: %d\n", i);

        // 2. Create a line, by specifying generated point and center of sphere/volume
        line l = {center, fibonacci_points[i]};

        std::vector<unsigned short> chord_lengths;
        unsigned long* result = chords_in_any_direction(vol, l, chords_through, chord_lengths);

        loose_chord_lengths[i] = chord_lengths;

        merge_chords(all_chords, result, num_elements, mode);
        free(result);

    }

    printf("Compressing image...");
    unsigned short* compressed = compress_results(all_chords, num_elements, number_of_angles, mode);
    printf("Done.\n");

    printf("Compressing chord lengths...");
    std::vector<unsigned short> chord_lengths;
    for (int i = 0; i < number_of_angles; i++) {
        chord_lengths.insert(chord_lengths.end(), loose_chord_lengths[i].begin(), loose_chord_lengths[i].end());
    }

    printf("Saving to file\n");
    // save to file final result
    std::string save_path = std::string(file_path) + "_" + std::to_string(number_of_angles) + "_" + mode;
    if (chords_through == 1) {
        save_path += "_inverse";
    }

    printf("number of chord lengths: %lu\n", chord_lengths.size());

    std::string image_save_path = save_path + ".raw";
    std::string lengths_save_path = save_path + ".dat";
    save_to_file(lengths_save_path, chord_lengths.data(), chord_lengths.size());
    save_to_file(image_save_path, compressed, num_elements);

    // Free the memory
    free(all_chords);
    free(compressed);

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
unsigned long* chords_in_any_direction(volume vol, line l, bool chords_through, std::vector<unsigned short>& chord_lengths) {
    unsigned long* chords = (unsigned long*) malloc((unsigned int) vol.nx * vol.ny * vol.nz * sizeof(unsigned long));

    // Check if memory allocation was successful
    if (chords == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return NULL;
    }

    // printf("Allocated memory\n");

    parametric_line pl = to_parametric_line(l);

    // printf("Converted to parametric line\n");

    unsigned int grid_size = vol.nx;
    parametric_line* parallel_lines = get_parallel_lines(pl, grid_size);

    //
    unsigned long num_paralell_lines = grid_size * grid_size * 3 * 4;
    for (unsigned int k = 0; k < num_paralell_lines; k++) {

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

        // DEBUG: Can't we stop this from happening
        if (stop.x >= vol.nx || stop.y >= vol.ny || stop.z >= vol.nz || start.x >= vol.nx || start.y >= vol.ny || start.z >= vol.nz) {
            //printf("Tried getting chords along line: %d, %d, %d, %d, %d, %d\n", start.x, start.y, start.z, stop.x, stop.y, stop.z);
            continue;
        }

        chords_along_line(vol, chords, start, stop, chords_through, chord_lengths);
    }

    free(parallel_lines);

    return chords;
}


// adjusted Bresenham's line algorithm to generate chords along a line and measure their length
void chords_along_line(volume vol, unsigned long* chords, point line_start, point line_end, bool chords_through, std::vector<unsigned short>& chord_lengths) {
    // FOR DEBUG ONLY
    // if (line_start.y % 2 == 1)
    //    return;

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
    point* start_point = 0;

    while (1) {  /* loop */
        unsigned long index = z0 * vol.ny * vol.nx + y0 * vol.nx + x0;
        unsigned short value = vol.data[index];

        // TODO: Chords through modus
        if (value == chords_through) {
                if (start_point == 0) {
                    start_point = (point*) malloc(sizeof(point));
                    start_point->x = static_cast<unsigned short>(x0);
                    start_point->y = static_cast<unsigned short>(y0);
                    start_point->z = static_cast<unsigned short>(z0);
                }
                chord_pixel* new_chord_pixel = (chord_pixel*) malloc(sizeof(chord_pixel));
                new_chord_pixel->previous_pixel = last_chord_pixel;
                new_chord_pixel->pixel_index = index;
                last_chord_pixel = new_chord_pixel;
        } else {
            if (last_chord_pixel != 0) {
                point end_point = {static_cast<unsigned short>(x0),static_cast<unsigned short>(y0),static_cast<unsigned short>(z0)};

                double distance = calc_distance(*start_point, end_point);
                unsigned short rounded_distance = round(distance);

                free(start_point);
                start_point = 0;

                set_chord(chords, *last_chord_pixel, distance);

                // add chord to chord lengths vector
                chord_lengths.push_back(distance);

                // Free memory
                while (1) {
                    chord_pixel* to_free = last_chord_pixel;
                    if (last_chord_pixel->previous_pixel != 0) {
                        last_chord_pixel = last_chord_pixel->previous_pixel;
                        free(to_free);
                    } else {
                        free(to_free);
                        last_chord_pixel = 0;
                        break;
                    }
                }
            }
        }

        if (i-- == 0) break;
        x1 -= dx; if (x1 < 0) { x1 += dm; x0 += sx; }
        y1 -= dy; if (y1 < 0) { y1 += dm; y0 += sy; }
        z1 -= dz; if (z1 < 0) { z1 += dm; z0 += sz; }
    }

    if (last_chord_pixel != 0) {
        point end_point = {static_cast<unsigned short>(x0),static_cast<unsigned short>(y0),static_cast<unsigned short>(z0)};
        double distance = calc_distance(*start_point, end_point);
		free(start_point);
        unsigned short rounded_distance = round(distance);
        set_chord(chords, *last_chord_pixel, rounded_distance);

		// Free memory

		while (1) {
            chord_pixel* to_free = last_chord_pixel;
            if (last_chord_pixel->previous_pixel != 0) {
                last_chord_pixel = last_chord_pixel->previous_pixel;
                free(to_free);
            } else {
                free(to_free);
                last_chord_pixel = 0;
                break;
            }
        }
    }

}

// use a stack to store chord pixels and set them efficiently
void set_chord(unsigned long* chords, chord_pixel last_pixel, unsigned int chord_length) {
    // Free up memory on the go, except for the last pixel, that will be freed by the caller
    while (1) {
        chords[last_pixel.pixel_index] = chord_length;
        if (last_pixel.previous_pixel != 0)
            last_pixel = *last_pixel.previous_pixel;
        else break;
    }
}

template<typename T>
void update_maximum(std::atomic<T>& maximum_value, T const& value) noexcept {
    T prev_value = maximum_value.load();
    while (prev_value < value && !maximum_value.compare_exchange_weak(prev_value, value)) {
        // If the compare_exchange_weak fails, prev_value is updated with the current value of maximum_value
    }
}

template<typename T>
void update_minimum(std::atomic<T>& minimum_value, T const& value) noexcept {
    T prev_value = minimum_value.load();
    while (prev_value > value && !minimum_value.compare_exchange_weak(prev_value, value)) {
        // If the compare_exchange_weak fails, prev_value is updated with the current value of minimum_value
    }
}

void merge_chords(unsigned long* chords, unsigned long* to_merge, unsigned long num_elements, std::string mode) {
    if (mode == "min") {
        for (int i = 0; i < num_elements; i++) {
            std::atomic<unsigned long> atomic_chord(chords[i]);
            update_minimum(atomic_chord, to_merge[i]);
            chords[i] = atomic_chord.load();
        }
    } else if (mode == "min") {
        for (int i = 0; i < num_elements; i++) {
            std::atomic<unsigned long> atomic_chord(chords[i]);
            update_maximum(atomic_chord, to_merge[i]);
            chords[i] = atomic_chord.load();
        }
    } else {
        for (int i = 0; i < num_elements; i++) {
            #pragma omp atomic
            chords[i] = chords[i] + to_merge[i];
        }
    }
}


// This function "compresses" the results => we work with longs internally, to be able to accomodate high accumulations in average mode
// In average mode the mean has to be taken here, in all other modes, we just cast and copy
unsigned short* compress_results(unsigned long* chords, unsigned long num_elements, unsigned short num_of_angles, std::string mode) {
    unsigned short* compressed = (unsigned short*) malloc((unsigned long) num_elements * sizeof(unsigned short));
    if (mode == "avg") {
        for (int i = 0; i < num_elements; i++) {
                compressed[i] = static_cast<unsigned short>(chords[i] / num_of_angles);
        }
    } else {
        for (int i = 0; i < num_elements; i++) {
            compressed[i] = static_cast<unsigned short>(chords[i]);
        }
    }
    return compressed;
}