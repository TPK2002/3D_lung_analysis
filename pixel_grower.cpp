#include <cfloat>
#include <cstdio>
#include <omp.h>
#include <sys/fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <cstring>
#include <string>
#include <limits.h>

#include "common.h"

#define MODE_MIN "min"
#define MODE_MAX "max"
#define MODE_AVERAGE "avg"

unsigned short get_pixel_distance(point current_point, volume vol, unsigned short number_of_angles, std::string mode, bool chords_through);
unsigned short get_line_length(volume vol, point line_start, point line_end, bool chords_through);
int save_to_file(std::string save_path, unsigned short* result, unsigned int num_elements);

int main(int argc, char *argv[]){

    // Get from CLI. nx, ny, nz, number_of_angles, mode
    std::string file_path;
    unsigned short nx, ny, nz, number_of_angles;
    std::string mode = "avg"; // default is average
    bool chords_through = 0; // should the processing happen on background: 0 or tissue: 1

    if (argc < 6) {
        printf("Pixel Grower has to be called with arguments: <file_path> <nx> <ny> <nz> <num_angles> (<mode: avg|min|max> < --inverse>)");
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

    printf("Pixel Grower called for: %s, %hux%hux%hu.\n", file_path.c_str(), nx, ny, nz);
    printf("Doing %hu angles per pixel with mode: %s \n",number_of_angles, mode.c_str());
    if (chords_through == 1) {
        printf("INVERSE MODE ACTIVE \n");
    }

    // Calculate the total number of elements
    unsigned int num_elements = nx * ny * nz;
    unsigned int file_size = num_elements;

    // Open the binary file
    int fd = open(file_path.c_str(), O_RDONLY);
    if (fd == -1) {
        perror("Error opening file");
        return 1;
    }

    // Get file properties, to compare to theoretical size
   	struct stat file_properties;
	if(fstat(fd, &file_properties) < 0)
	{
		perror("Error getting file properties");
		return 1;
	}

    // Compare theoretical and practical file size and quit if they don't align
	unsigned int N_bytes = file_properties.st_size;
	if ((unsigned int) N_bytes != (unsigned int) file_size)
    {
        printf("File size: %u\n", (unsigned int) N_bytes);
        printf("Expected size: %u\n", (unsigned int) file_size);
        printf("File size does not match the expected size\n");
        return 1;
    }

    // Map the file into the processes virtual memory address space
    unsigned char *data = (unsigned char *)mmap(NULL, (unsigned int) file_size, PROT_READ, MAP_SHARED, fd, 0);
    if (data == MAP_FAILED) {
        perror("Error mmapping the file");
        close(fd);
        return 1;
    }


    // All is prepared, we can begin with the actual processing
    printf("Beginning processing\n");

    // Create volume structure, for easier access to parameters
    volume vol = {nx,ny,nz,data};

    // Reserve Memory for results
    unsigned short* result = (unsigned short*) malloc((unsigned int) vol.nx * vol.ny * vol.nz * sizeof(unsigned short));

    // This little preprocessor directive does all the magic needed for parallel processing
    # pragma omp parallel for

    // Do for every slice
    for (unsigned short z = 0; z < vol.nz; z++) {

        // Do for every row
        for (unsigned short y = 0; y < vol.ny; y++) {

            // Do for every pixel
            for (unsigned short x = 0; x < vol.nx; x++) {

                // Calculate pixel index of currently looked at point
                point current_point = {x,y,z};
                unsigned long pixel_index = z * vol.ny * vol.nx + y * vol.nx + x;

                // We only do calculations for pixels, that interest us, either background or tissue
                if (vol.data[pixel_index] == chords_through)
                    result[pixel_index] = get_pixel_distance(current_point, vol, number_of_angles, mode, chords_through);
            }
        }
    }

    std::string save_path = std::string(file_path) + "_" + std::to_string(nx) + "x" + std::to_string(ny) + "x" + std::to_string(nz) + "_" + std::to_string(number_of_angles) + "_" + mode;
    if (chords_through == 1) {
        save_path += "_inverse";
    }
    save_path += ".raw";

    // save to file final result
    save_to_file(save_path, result, num_elements);

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

unsigned short get_pixel_distance(point current_point, volume vol, unsigned short number_of_angles, std::string mode, bool chords_through){
    // 1. Calculate fibonacci Sphere Points
    // CAN WE OPTIMIZE, BY NOT RECALCULATING THE SPHERE WITH EVERY PIXEL ?
    point fibonacci_points[number_of_angles];
    unsigned short fibonacci_radius = min(vol.nx, vol.ny, vol.nz) / 2;
    fibonacci_sphere(fibonacci_points, number_of_angles, current_point, fibonacci_radius);

    unsigned short min_distance = USHRT_MAX;
    unsigned short max_distance = 0;
    unsigned short cumulated_distance = 0;
    for (int i = 0;i < number_of_angles; i++) {

        // Get max line bound
        double t_min, t_max;
        parametric_line pl = to_parametric_line({current_point, fibonacci_points[i]});
        get_line_bounds(pl, vol, &t_min, &t_max);

        // TODO: Make sure this cast does not produce values > volume size
        point end_point = {static_cast<unsigned short>(pl.x0 + pl.direction.a * t_max), static_cast<unsigned short>(pl.y0 + pl.direction.b * t_max), static_cast<unsigned short>(pl.z0 + pl.direction.c * t_max)};

        // Get line length until we hit a non-zero voxel or the end of the line
        unsigned short line_length = get_line_length(vol, current_point, end_point, chords_through);
        cumulated_distance += line_length;
        if (line_length < min_distance) min_distance = line_length;
        if (line_length > max_distance) max_distance = line_length;
    }

    // Depending on the mode combine the values:
    if (mode == MODE_MIN) return min_distance;
    if (mode == MODE_MAX) return max_distance;
    if (mode == MODE_AVERAGE) return cumulated_distance / number_of_angles;

    fprintf(stderr, "Unknown mode\n");
    exit(EXIT_FAILURE);
}

// adjusted Bresenham's line algorithm to walk along a line until we hit a non-zero/non-one voxel
unsigned short get_line_length(volume vol, point line_start, point line_end, bool chords_through) {
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

    unsigned int chord_length = 0;

    while (1) {  /* loop */
        unsigned int pixel_index = z0 * vol.nx * vol.ny + y0 * vol.nx + x0;
        unsigned short value = vol.data[pixel_index];

        // Move on the line until we hit a non-zero voxel
        if (value != chords_through) {
            return chord_length;
        }
        chord_length++;

        if (i-- == 0) break;
        x1 -= dx; if (x1 < 0) { x1 += dm; x0 += sx; }
        y1 -= dy; if (y1 < 0) { y1 += dm; y0 += sy; }
        z1 -= dz; if (z1 < 0) { z1 += dm; z0 += sz; }
    }

    return chord_length;
}

int save_to_file(std::string save_path, unsigned short* result, unsigned int num_elements) {
    // save to file
    FILE *fp = fopen(save_path.c_str(), "wb");
    if (fp == NULL) {
        perror("Error opening file");
        return 1;
    }
    fwrite(result, sizeof(unsigned short), num_elements, fp);
    fclose(fp);

    return 0;
}
