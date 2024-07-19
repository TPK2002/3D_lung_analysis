#include <cfloat>
#include <cmath>
#include <cstdio>

struct volume {
    unsigned short nx;
    unsigned short ny;
    unsigned short nz;
    unsigned char* data;
};

struct point {
    unsigned short x;
    unsigned short y;
    unsigned short z;
};

struct vector {
    double a;
    double b;
    double c;
};

struct line {
    point start;
    point end;
};

struct parametric_line {
    double x0;
    double y0;
    double z0;
    vector direction;
};

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

vector get_perpendicular_vector(vector l) {
    vector perp;
    perp.a = -l.b;
    perp.b = l.a;
    perp.c = 0;
    return normalize(perp);
}

vector get_crossprod_perpendicular(vector v1, vector v2) {
    vector perpendicular;
    perpendicular.a = v1.b * v2.c - v1.c * v2.b;
    perpendicular.b = v1.c * v2.a - v1.a * v2.c;
    perpendicular.c = v1.a * v2.b - v1.b * v2.a;
    return normalize(perpendicular);
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

void get_line_bounds(parametric_line line, volume vol, double *t_min, double *t_max) {
    double x_min = 0, y_min = 0, z_min = 0;
    double x_max = vol.nx - 1, y_max = vol.ny - 1, z_max = vol.nz - 1;

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

parametric_line* get_parallel_lines(parametric_line l, unsigned short grid_size) {
    parametric_line* parallel_lines = (parametric_line*) malloc(sizeof(parametric_line) * grid_size * grid_size * 3 * 4);
    if (parallel_lines == NULL) {
        perror("Failed to allocate memory for parallel lines\n");
        return NULL;
    }

    vector perp1 = get_perpendicular_vector(l.direction);
    vector perp2 = get_crossprod_perpendicular(l.direction, perp1);

    // TODO: explain formula:
    // https://www.mathe-lexikon.at/geometrie/geometrische-koerper/wuerfel/raumdiagonale.html#google_vignette
    int a = sqrt(3) * grid_size;
    for (int i = -a; i < a; i++) {
        for (int j = -a; j < a; j++) {
            // We intentionally include the center line
            vector shift = {i/2.0 * perp1.a + j/2.0 * perp2.a, i/2.0 * perp1.b + j/2.0 * perp2.b, i/2.0 * perp1.c + j/2.0 * perp2.c};
            parametric_line new_line;
            new_line.x0 = l.x0 + shift.a;
            new_line.y0 = l.y0 + shift.b;
            new_line.z0 = l.z0 + shift.c;
            new_line.direction = l.direction;

            // TODO: Explain index calculation
            unsigned long index = (i + a) * (a*2) + (j + a);
            parallel_lines[index] = new_line;
        }
    }

    // DEBUG
    // // Print the parallel lines
    // for (int i = 0; i < grid_size * grid_size * 2 * 2; i++) {
    //     printf("Line %d: x0: %f, y0: %f, z0: %f, a: %f, b: %f, c: %f\n", i, parallel_lines[i].x0, parallel_lines[i].y0, parallel_lines[i].z0, parallel_lines[i].direction.a, parallel_lines[i].direction.b, parallel_lines[i].direction.c);
    // }

    return parallel_lines;
}

// Hypot is highly optimized to e.g. avoid overflow
double calc_distance(point start, point end) {
    return hypot(hypot(start.x-end.x,start.y-end.y),start.z-end.z);
}

void fibonacci_hemisphere(point* fibonacci_points, unsigned int number_of_angles, point center, double radius) {
    double golden_ratio = (1.0 + sqrt(5.0)) / 2.0;

    double n = number_of_angles;

    for (int i = 0; i < n; i++) {
        double theta = 2 * M_PI * i / golden_ratio;

        // i + 0.5 verbessert die Genauigkeit wesentlich!
        // ERKLÄREN
        double p = 2 * n + 1;
        double phi = acos((2 * i) / p);

        double x = cos(theta) * sin(phi);
        double y = sin(theta) * sin(phi);
        double z = cos(phi);

        fibonacci_points[i].x = center.x + radius * x;
        fibonacci_points[i].y = center.y + radius * y;
        fibonacci_points[i].z = center.z + radius * z;
    }
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

template<typename T>
int save_to_file(std::string save_path, T* result, unsigned int num_elements) {
    // save to file
    FILE *fp = fopen(save_path.c_str(), "wb");
    if (fp == NULL) {
        perror("Error opening file");
        return 1;
    }
    fwrite(result, sizeof(T), num_elements, fp);
    fclose(fp);

    return 0;
}