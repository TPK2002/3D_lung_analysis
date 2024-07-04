#include <cfloat>
#include <cstdio>
#include <omp.h>
#include <cmath>
#include <sys/fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

// either 0 or 1, depending on what we to end a chord with
// We will also only work on pixels, that are not equal to CHORDS_UNTIL
bool CHORDS_UNTIL = 1;

struct point {
    unsigned short x;
    unsigned short y;
    unsigned short z;
};

struct volume {
    unsigned short width;
    unsigned short height;
    unsigned short depth;
    unsigned char* data;
};

struct line {
    point start;
    point end;
};

struct vector {
    double a;
    double b;
    double c;
};

struct parametric_line {
    double x0;
    double y0;
    double z0;
    vector direction;
};

int max(int a, int b, int c) {
    if (a > b && a > c) {
        return a;
    } else if (b > a && b > c) {
        return b;
    } else {
        return c;
    }
}

vector normalize(vector v) {
    double length = sqrt(v.a * v.a + v.b * v.b + v.c * v.c);
    vector normalized;
    normalized.a = v.a / length;
    normalized.b = v.b / length;
    normalized.c = v.c / length;
    return normalized;
}

// adjusted Bresenham's line algorithm to walk along a line until we hit a non-zero/non-one voxel
unsigned short get_line_length(volume vol, point line_start, point line_end) {
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
        unsigned int index = z0 * vol.width * vol.height + x0 * vol.height + y0;
        unsigned short value = vol.data[index];

        // Move on the line until we hit a non-zero voxel
        if (value == CHORDS_UNTIL) {
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

void get_line_bounds(parametric_line line, volume vol, double *t_min, double *t_max) {
    double x_min = 0, y_min = 0, z_min = 0;
    double x_max = vol.width - 1, y_max = vol.height - 1, z_max = vol.depth - 1;

    double t_x_min, t_x_max;
    double t_y_min, t_y_max;
    double t_z_min, t_z_max;

    // Avoid division by zero by handling a, b, or c being zero
    if (line.direction.a != 0) {
        t_x_min = (x_min - line.x0) / line.direction.a;
        t_x_max = (x_max - line.x0) / line.direction.a;
        if (t_x_min > t_x_max) {
            double temp = t_x_min;
            t_x_min = t_x_max;
            t_x_max = temp;
        }
    } else {
        t_x_min = -DBL_MAX;
        t_x_max = DBL_MAX;
    }

    if (line.direction.b != 0) {
        t_y_min = (y_min - line.y0) / line.direction.b;
        t_y_max = (y_max - line.y0) / line.direction.b;
        if (t_y_min > t_y_max) {
            double temp = t_y_min;
            t_y_min = t_y_max;
            t_y_max = temp;
        }
    } else {
        t_y_min = -DBL_MAX;
        t_y_max = DBL_MAX;
    }

    if (line.direction.c != 0) {
        t_z_min = (z_min - line.z0) / line.direction.c;
        t_z_max = (z_max - line.z0) / line.direction.c;
        if (t_z_min > t_z_max) {
            double temp = t_z_min;
            t_z_min = t_z_max;
            t_z_max = temp;
        }
    } else {
        t_z_min = -DBL_MAX;
        t_z_max = DBL_MAX;
    }

    // Find the intersection of all valid t ranges
    *t_min = fmax(fmax(t_x_min, t_y_min), t_z_min);
    *t_max = fmin(fmin(t_x_max, t_y_max), t_z_max);

    if (*t_min > *t_max) {
        // ------------ TODO: Investigate this, wasted computing power, no ? ------------
        //printf("No valid t range: t_min = %f, t_max = %f\n", *t_min, *t_max);
        *t_min = -1; // Or some indication of no valid range
        *t_max = -1;
    }
}

parametric_line to_parametric_line(line l) {
    parametric_line pl;
    pl.x0 = l.start.x;
    pl.y0 = l.start.y;
    pl.z0 = l.start.z;
    pl.direction.a = l.end.x - l.start.x;
    pl.direction.b = l.end.y - l.start.y;
    pl.direction.c = l.end.z - l.start.z;
    pl.direction = normalize(pl.direction);
    return pl;
}

void fibonacci_sphere(point* fibonacci_points, unsigned int number_of_angles, point center, double radius) {
    double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;

    for (int i = 0; i < number_of_angles; i++) {
        double theta = 2.0 * M_PI * i / golden_ratio;
        double phi = acos(1 - 2 * i / (double) number_of_angles); // eventuell i + 0,5

        double x = cos(theta) * sin(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(phi);

        fibonacci_points[i].x = center.x + radius * x;
        fibonacci_points[i].y = center.y + radius * y;
        fibonacci_points[i].z = center.z + radius * z;
    }
}

int min(int a,int b,int c) {
    if (a < b) {
        if (a < c) {
            return a;
        } else {
            return c;
        }
    } else {
        if (b < c) {
            return b;
        } else {
            return c;
        }
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

unsigned short get_pixel_distance(point current_point, volume vol, unsigned short number_of_angles){

    // 1. Calculate fibonacci Sphere Points
    point fibonacci_points[number_of_angles];
    unsigned short fibonacci_radius = min(vol.width, vol.height, vol.depth) / 2;
    fibonacci_sphere(fibonacci_points, number_of_angles, current_point, fibonacci_radius);

    unsigned long cumulated_distance = 0;

    # pragma omp parallel for
    for (int i = 0;i < number_of_angles; i++) {

        // Get max line bound
        double t_min, t_max;
        parametric_line pl = to_parametric_line({current_point, fibonacci_points[i]});
        get_line_bounds(pl, vol, &t_min, &t_max);

        // TODO: Make sure this cast does not produce values > volume size
        point end_point = {static_cast<unsigned short>(pl.x0 + pl.direction.a * t_max), static_cast<unsigned short>(pl.y0 + pl.direction.b * t_max), static_cast<unsigned short>(pl.z0 + pl.direction.c * t_max)};

        // Get line length until we hit a non-zero voxel or the end of the line
        cumulated_distance += get_line_length(vol, current_point, end_point);
    }

    // return mean value
    return cumulated_distance / number_of_angles;
}


int main(){
    // Dimensions of the 3D array
    int width = 500;  // replace with actual width
    int height = 500; // replace with actual height
    int depth = 500;  // replace with actual depth
    unsigned short number_of_angles = 1000;

    // Calculate the total number of elements
    unsigned int num_elements = width * height * depth;
    unsigned int file_size = num_elements;

    // Open the binary file
    int fd = open("./lung424_small_mask.raw", O_RDONLY);
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

    unsigned long* result = (unsigned long*) malloc((unsigned int) vol.width * vol.height * vol.depth * sizeof(unsigned long));

    // For Every Pixel do:
    for (unsigned short i = 0; i < vol.depth; i++) {
        printf("Working on slice no.: %d\n", i);
        for (unsigned short j = 0; j < vol.width; j++) {
            for (unsigned short k = 0; k < vol.height; k++) {
                point current_point = {j,k,i};
                unsigned long pixel_index = i * vol.width * vol.height + j * vol.height + k;
                result[pixel_index] = get_pixel_distance(current_point, vol, number_of_angles);
            }
            if (j % 20 == 0) {
                save_to_file(result, num_elements);
            }
        }
    }

    // save to file final result
    save_to_file(result, num_elements);

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
