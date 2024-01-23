// Author: APD team, except where source was noted

#include "helpers.h"
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <pthread.h>

#define CONTOUR_CONFIG_COUNT    16
#define FILENAME_MAX_SIZE       50
#define STEP                    8
#define SIGMA                   200
#define RESCALE_X               2048
#define RESCALE_Y               2048

#define CLAMP(v, min, max) if(v < min) { v = min; } else if(v > max) { v = max; }
#define MIN(a, b) (a < b ? a : b)

typedef struct {
    int thread_id;
    int nr_threads;
    unsigned char **grid;
    int step_x;
    int step_y;
    ppm_image *image;
    ppm_image *old_image;
    ppm_image **contour_map;
    pthread_barrier_t *barrier;
} thread_data;

// Creates a map between the binary configuration (e.g. 0110_2) and the corresponding pixels
// that need to be set on the output image. An array is used for this map since the keys are
// binary numbers in 0-15. Contour images are located in the './contours' directory.
ppm_image **init_contour_map() {
    ppm_image **map = (ppm_image **)malloc(CONTOUR_CONFIG_COUNT * sizeof(ppm_image *));
    if (!map) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        char filename[FILENAME_MAX_SIZE];
        sprintf(filename, "./contours/%d.ppm", i);
        map[i] = read_ppm(filename);
    }

    return map;
}

// Updates a particular section of an image with the corresponding contour pixels.
// Used to create the complete contour image.
void update_image(ppm_image *image, ppm_image *contour, int x, int y) {
    for (int i = 0; i < contour->x; i++) {
        for (int j = 0; j < contour->y; j++) {
            int contour_pixel_index = contour->x * i + j;
            int image_pixel_index = (x + i) * image->y + y + j;

            image->data[image_pixel_index].red = contour->data[contour_pixel_index].red;
            image->data[image_pixel_index].green = contour->data[contour_pixel_index].green;
            image->data[image_pixel_index].blue = contour->data[contour_pixel_index].blue;
        }
    }
}

void free_resources(ppm_image *image, ppm_image **contour_map, unsigned char **grid, int step_x) {
    for (int i = 0; i < CONTOUR_CONFIG_COUNT; i++) {
        if (contour_map[i]->data != NULL)
            free(contour_map[i]->data);
        
        if (contour_map[i] != NULL)
            free(contour_map[i]);
    }
    if (contour_map != NULL)
        free(contour_map);

    for (int i = 0; i <= image->x / step_x; i++) {
        if (grid[i] != NULL)
            free(grid[i]);
    }
    if (grid != NULL)
        free(grid);

    if (image->data != NULL)
        free(image->data);

    if (image != NULL)
        free(image);
}


void *thread_function(void *arg) {
    thread_data *data = (thread_data *)arg;

    // rescale image if needed
    if (data->old_image->x > RESCALE_X || data->old_image->y > RESCALE_Y) {
        uint8_t sample[3];

        int start = data->thread_id * data->image->x * data->image->y / data->nr_threads;
        int end = MIN((data->thread_id + 1) * data->image->x * data->image->y / data->nr_threads, data->image->x * data->image->y);

        // use bicubic interpolation for scaling
        for (int i = start; i < end; i++) {
            int x = i / data->image->y;
            int y = i % data->image->y;

            float u = (float)x / (float)(data->image->x - 1);
            float v = (float)y / (float)(data->image->y - 1);
            sample_bicubic(data->old_image, u, v, sample);

            data->image->data[i].red = sample[0];
            data->image->data[i].green = sample[1];
            data->image->data[i].blue = sample[2];
        }

    }

    //bariera
    pthread_barrier_wait(data->barrier);

    // Corresponds to step 1 of the marching squares algorithm, which focuses on sampling the image.
    // Builds a p x q grid of points with values which can be either 0 or 1, depending on how the
    // pixel values compare to the `sigma` reference value. The points are taken at equal distances
    // in the original image, based on the `step_x` and `step_y` arguments.
    int p = data->image->x / data->step_x;
    int q = data->image->y / data->step_y;

    // set starting and ending points for each thread
    int start = data->thread_id * p * q / data->nr_threads;
    int end = MIN((data->thread_id + 1) * p * q / data->nr_threads, p * q);

    for (int i = start; i < end; i++) {
        int x = i / q;
        int y = i % q;

        ppm_pixel curr_pixel = data->image->data[x * data->step_x * data->image->y + y * data->step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            data->grid[x][y] = 0;
        } else {
            data->grid[x][y] = 1;
        }
    }

    data->grid[p][q] = 0;

    // last sample points have no neighbors below / to the right, so we use pixels on the
    // last row / column of the input image for them
    start = data->thread_id * p / data->nr_threads;
    end = MIN((data->thread_id + 1) * p / data->nr_threads, p);
    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = data->image->data[i * data->step_x * data->image->y + data->image->x - 1];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            data->grid[i][q] = 0;
        } else {
            data->grid[i][q] = 1;
        }
    }

    // set starting and ending points for each thread
    start = data->thread_id * q / data->nr_threads;
    end = MIN((data->thread_id + 1) * q / data->nr_threads, q);
    for (int i = start; i < end; i++) {
        ppm_pixel curr_pixel = data->image->data[(data->image->x - 1) * data->image->y + i * data->step_y];

        unsigned char curr_color = (curr_pixel.red + curr_pixel.green + curr_pixel.blue) / 3;

        if (curr_color > SIGMA) {
            data->grid[p][i] = 0;
        } else {
            data->grid[p][i] = 1;
        }
    }

    // wait for all threads to finish sampling
    pthread_barrier_wait(data->barrier);


    // set starting and ending points for each thread for marching
    start = data->thread_id * p * q / data->nr_threads;
    end = MIN((data->thread_id + 1) * p * q / data->nr_threads, p * q);

    // Corresponds to step 2 of the marching squares algorithm, which focuses on identifying the
    // type of contour which corresponds to each subgrid. It determines the binary value of each
    // sample fragment of the original image and replaces the pixels in the original image with
    // the pixels of the corresponding contour image accordingly.
    for (int i = start; i < end; i++) {
        int x = i / q;
        int y = i % q;

        // update image
        unsigned char k = 8 * data->grid[x][y] + 4 * data->grid[x][y + 1] + 2 * data->grid[x + 1][y + 1] + 1 * data->grid[x + 1][y];
        update_image(data->image, data->contour_map[k], x * data->step_x, y * data->step_y);
    }

    pthread_exit(NULL);
}

int main(int argc, char *argv[]) {
    if (argc < 4) {
        fprintf(stderr, "Usage: ./tema1 <in_file> <out_file> <P>\n");
        return 1;
    }

    int num_threads = atoi(argv[3]);
    ppm_image *image = read_ppm(argv[1]);

    ppm_image **contour_map = init_contour_map();
    ppm_image *scaled_image;

    if (image->x <= RESCALE_X && image->y <= RESCALE_Y) {
        scaled_image = image;
    } else {
        // allocate memory for scaled image
        scaled_image = (ppm_image *)malloc(sizeof(ppm_image));

        if (!scaled_image) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }

        scaled_image->x = RESCALE_X;
        scaled_image->y = RESCALE_Y;

        scaled_image->data = (ppm_pixel*)malloc(scaled_image->x * scaled_image->y * sizeof(ppm_pixel));
        if (!scaled_image->data) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    int step_x = STEP;
    int step_y = STEP;
    // allocate memory for grid
    unsigned char **grid = (unsigned char **)malloc((scaled_image->x / step_x + 1) * sizeof(unsigned char*));
    if (!grid) {
        fprintf(stderr, "Unable to allocate memory\n");
        exit(1);
    }

    for (int i = 0; i <= scaled_image->x / step_x; i++) {
        grid[i] = (unsigned char *)malloc((scaled_image->y / step_y + 1) * sizeof(unsigned char));
        if (!grid[i]) {
            fprintf(stderr, "Unable to allocate memory\n");
            exit(1);
        }
    }

    // create barrier
    pthread_barrier_t *barrier = (pthread_barrier_t *)malloc(sizeof(pthread_barrier_t));
    pthread_barrier_init(barrier, NULL, num_threads);

    // initialize thread data vector
    thread_data data[num_threads];
    for (int i = 0; i < num_threads; i++) {
        data[i].thread_id = i;
        data[i].nr_threads = num_threads;
        data[i].grid = grid;
        data[i].step_x = step_x;
        data[i].step_y = step_y;
        data[i].image = scaled_image;
        data[i].old_image = image;
        data[i].contour_map = contour_map;
        data[i].barrier = barrier;
    }

    // Create threads
    pthread_t tid[num_threads];

    for (int i = 0; i < num_threads; i++) {
        printf("Creating thread %d\n", i);
        pthread_create(&(tid[i]), NULL, thread_function, &(data[i]));
    }

    for (int i = 0; i < atoi(argv[3]); i++) {
        pthread_join(tid[i], NULL);
    }

    pthread_barrier_destroy(data[0].barrier);

    // 4. Write output
    write_ppm(data[0].image, argv[2]);

    return 0;
}
