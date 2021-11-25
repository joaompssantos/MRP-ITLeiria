#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mrp.h"

//Reference pixel position separated by their distances {y, x}
const POINT dyx[] = {
        /* 1 --> 2 */
        { 0,-1}, {-1, 0},
        /* 2 --> 6 */
        { 0,-2}, {-1,-1}, {-2, 0}, {-1, 1},
        /* 3 --> 12 */
        { 0,-3}, {-1,-2}, {-2,-1}, {-3, 0}, {-2, 1}, {-1, 2},
        /* 4 --> 20 */
        { 0,-4}, {-1,-3}, {-2,-2}, {-3,-1}, {-4, 0}, {-3, 1}, {-2, 2}, {-1, 3},
        /* 5 --> 30 */
        { 0,-5}, {-1,-4}, {-2,-3}, {-3,-2}, {-4,-1}, {-5, 0}, {-4, 1}, {-3, 2}, {-2, 3}, {-1, 4},
        /* 6 --> 42 */
        { 0,-6}, {-1,-5}, {-2,-4}, {-3,-3}, {-4,-2}, {-5,-1}, {-6, 0}, {-5, 1}, {-4, 2}, {-3, 3}, {-2, 4}, {-1, 5},
        /* 7 --> 56 */
        { 0,-7}, {-1,-6}, {-2,-5}, {-3,-4}, {-4,-3}, {-5,-2}, {-6,-1}, {-7, 0}, {-6, 1}, {-5, 2}, {-4, 3}, {-3, 4}, {-2, 5}, {-1, 6},
        /* 8 --> 72 */
        { 0,-8}, {-1,-7}, {-2,-6}, {-3,-5}, {-4,-4}, {-5,-3}, {-6,-2}, {-7,-1}, {-8, 0}, {-7, 1}, {-6, 2}, {-5, 3}, {-4, 4}, {-3, 5}, {-2, 6}, {-1, 7},
};

//Reference pixel position separated by their distances
const POINT idyx[] = {
        /* 0 --> 1 */
        //{ 0, 0},
        /* 1 --> 5 */
        { 0,-1}, {-1, 0}, { 0, 1}, { 1, 0},
        /* 2 --> 13 */
        { 0,-2}, {-1,-1}, {-2, 0}, {-1, 1}, { 0, 2}, { 1, 1}, { 2, 0}, { 1,-1},
        /* 3 --> 25 */
        { 0,-3}, {-1,-2}, {-2,-1}, {-3, 0}, {-2, 1}, {-1, 2}, { 0, 3}, { 1, 2}, { 2, 1}, { 3, 0}, { 2,-1}, { 1,-2},
        /* 4 --> 41 */
        { 0,-4}, {-1,-3}, {-2,-2}, {-3,-1}, {-4, 0}, {-3, 1}, {-2, 2}, {-1, 3}, { 0, 4}, { 1, 3}, { 2, 2}, { 3, 1}, { 4, 0}, { 3,-1}, { 2,-2}, { 1,-3},
};

const int unsigned mask[9] = {0x00FF,
                              0x01FF,
                              0x03FF,
                              0x07FF,
                              0x0FFF,
                              0x1FFF,
                              0x3FFF,
                              0x7FFF,
                              0xFFFF
};

double sigma_h[] = {0.85, 1.15, 1.50, 1.90, 2.55, 3.30, 4.25, 5.60, 7.15, 9.20, 12.05, 15.35, 19.95, 25.85, 32.95, 44.05};
double sigma_a[] = {0.15, 0.26, 0.38, 0.57, 0.83, 1.18, 1.65, 2.31, 3.22, 4.47, 6.19, 8.55, 11.80, 16.27, 22.42, 30.89};
double qtree_prob[] = {0.05, 0.2, 0.35, 0.5, 0.65, 0.8, 0.95};

/* Alternative version for 'free()' */
void safefree(void **pp) {
    /* in debug mode, abort if pp is NULL */
    assert(pp);
    if (pp != NULL) {               /* safety check */
        free(*pp);                  /* deallocate chunk, note that free(NULL) is valid */
        *pp = NULL;                 /* reset original pointer */
    }
}

void safefree_lf4d(LF4D **pp) {
    /* in debug mode, abort if pp is NULL */
    assert(pp);
    if (pp != NULL) {               /* safety check */
        free((*pp)->val);
        free(*pp);                  /* deallocate chunk, note that free(NULL) is valid */
        *pp = NULL;                 /* reset original pointer */
    }
}

/*------------------------------- fileopen --------------------------*
 |  Function fileopen
 |
 |  Purpose:  Opens a file, and verifies if it was open
 |
 |  Parameters:
 |      filename 	--> Name of the file (IN)
 |		mode		--> Opening mode (IN)
 |
 |  Returns:  fp* --> returns a file type structure
 *-------------------------------------------------------------------*/
FILE *fileopen(char *filename, char *mode) {
    FILE *fp;
    fp = fopen(filename, mode);
    if (fp == NULL) {
        fprintf(stderr, "Can\'t open %s!\n", filename);
        exit(1);
    }
    return (fp);
}

// Memory allocation
void *alloc_mem(size_t size) {
    void *ptr;

    // Check if memory was allocated
    if ((ptr = malloc(size)) == NULL) {
        fprintf(stderr, "Can\'t allocate memory (size = %d)!\n", (int)size);
        exit(1);
    }
    return (ptr);
}

// Image matrix allocation
void **alloc_2d_array(size_t height, size_t width, size_t size) {
    void **mat;
    char *ptr;
    size_t k;

    mat = (void **) alloc_mem(sizeof(void *) * height + height * width * size);
    ptr = (char *) (mat + height);

    for (k = 0; k < height; k++) {
        mat[k] =  ptr;
        ptr += width * size;
    }

    return (mat);
}

// Image matrix allocation
void ***alloc_3d_array(size_t frames, size_t height, size_t width, size_t size) {
    void ***mat;
    size_t i, j;

    mat = (void ***) alloc_mem(sizeof(void **) * frames + sizeof(void *) * frames * height + size * frames * height * width);

    void **hei = (void **) (mat + frames);
    char *wid = (char *) (hei + frames * height);

    for (i = 0; i < frames; i++) {
        mat[i] = &hei[i * height];
    }

    for (i = 0; i < frames; i++) {
        for (j = 0; j < height; j++) {
            hei[i * height + j] = &wid[i * height * width * size + j * width * size];
        }
    }

    return (mat);
}

// 4D matrix allocation
void ****alloc_4d_array(size_t v, size_t u, size_t t, size_t s, size_t size) {
    void ****mat;
    unsigned int i, j, k;

    // mat = (void ***) alloc_mem(sizeof(void ***) * frames + sizeof(void **) * frames * height + size * frames * height * width);
    mat = (void ****) alloc_mem(sizeof(void ***) * v + sizeof(void **) * v * u + sizeof(void *) * v * u * t + size * v * u * t * s);

    void ***vv = (void ***) (mat + v);
    void **vu = (void **) (vv + v * u);
    char *vut = (char *) (vu + v * u * t);

    for (i = 0; i < v; i++) {
        mat[i] = vv;
        vv += u;

        for (j = 0; j < u; j++) {
            mat[i][j] = vu;
            vu += t;

            for (k = 0; k < t; k++) {
                mat[i][j][k] = vut;
                vut += s * size;
            }
        }
    }

    return (mat);
}

// Function to alloc image
LF4D *alloc_lf4d( int v, int u, int t, int s, int maxval) {
    // Image struct vector
    LF4D *lf;

    // Image struct memory allocation
    lf = (LF4D *) alloc_mem(sizeof(LF4D));

    // Set image values
    lf->t = t;
    lf->s = s;
    lf->v = v;
    lf->u = u;

    lf->maxval = maxval;

    lf->val = (img_t ****) alloc_4d_array(lf->v, lf->u, lf->t, lf->s, sizeof(img_t));

    return (lf);
}

// Function to reverse the endianness of a unsigned short
unsigned short reverse_endianness (unsigned short s, int endianness) {
    unsigned char c1, c2;

    if (endianness == BIG_ENDIANNESS) {
        return s;
    }
    else {
        c1 = (unsigned char) (s & 255);
        c2 = (unsigned char) ((s >> 8) & 255);

        return (c1 << 8) + c2;
    }
}

// Write YUV image to file
void write_yuv(LF4D *lf, char *filename, int depth, int endianness, int format) {
    int v, u, t, s;
    unsigned short byte;
    FILE *fp;

    fp = fileopen(filename, "ab");

    int elements = lf->v * lf->u * lf->t * lf->s * (depth == 8 ? 1 : 2);
    unsigned char *stream_ptr, *stream = (unsigned char *) alloc_mem(elements * sizeof(unsigned char));
    stream_ptr = stream;

    switch (format) {
        case PVS:
            for (v = 0; v < lf->v; v++) {
                for (u = 0; u < lf->u; u++) {
                    for (t = 0; t < lf->t; t++) {
                        for (s = 0; s < lf->s; s++) {
                            if (depth == 8) {
                                *stream_ptr++ = lf->val[v][u][t][s];
                            }
                            else if (depth > 8) {
                                byte = reverse_endianness(lf->val[v][u][t][s], endianness);

                                *stream_ptr++ = (byte >> 8u) & 0x00FFu;
                                *stream_ptr++ = byte & 0x00FFu;
                            }
                        }
                    }
                }
            }

            break;

        case SAI:
        default:
            for (v = 0; v < lf->v; v++) {
                for (t = 0; t < lf->t; t++) {
                    for (u = 0; u < lf->u; u++) {
                        for (s = 0; s < lf->s; s++) {
                            if (depth == 8) {
                                *stream_ptr++ = lf->val[v][u][t][s];
                            }
                            else if (depth > 8) {
                                byte = reverse_endianness(lf->val[v][u][t][s], endianness);

                                *stream_ptr++ = (byte >> 8u) & 0x00FFu;
                                *stream_ptr++ = byte & 0x00FFu;
                            }
                        }
                    }
                }
            }

            break;
    }

    fwrite(stream, sizeof(unsigned char), elements, fp);

    free(stream);
    fclose(fp);
}

/*------------------------------- read_yuv --------------------------*
 |  Function read_yuv
 |
 |  Purpose:  Reads a yuv file to the program memory
 |
 |  Parameters:
 |      filename 	--> Name of the yuv file (IN)
 |		height		--> Height of the video (IN)
 |		width		--> Width of the video (IN)
 |		frame		--> Frame of the video to copy (IN)
 |		depth		--> Image dtnamic range (bpp) (IN)
 |		endianness	--> YUV endianness, for depth > 8 bpp (IN)
 |		format      --> Light field file format (IN)
 |
 |  Returns:  IMAGE* --> returns a video type structure
 *-------------------------------------------------------------------*/
LF4D *read_yuv(char *filename, int v, int u, int t, int s, int depth, int endianness, int format) {
    int i, j, k, l;
    unsigned int shift_first, shift_second, first, second;
    LF4D *lf;
    FILE *fp;


    unsigned long elements = v * u * t * s * (depth == 8 ? 1 : 2);
    unsigned char *stream_ptr, *stream = (unsigned char *) alloc_mem(elements * sizeof(unsigned char));

    //Open file
    fp = fileopen(filename, "rb");

    if (fread(stream, sizeof(unsigned char), elements, fp) != elements) {
        fprintf(stderr, "An error occurred while reading the image from file.\n");
        exit(EXIT_FAILURE);
    }
    stream_ptr = stream;

    if (endianness == LITTLE_ENDIANNESS) {
        shift_first = 0;
        shift_second = 8;
    }
    else {
        shift_first = 8;
        shift_second = 0;
    }

    // Image allocation
    lf = alloc_lf4d(v, u, t, s, (int) (pow(2, depth) - 1));

    switch (format) {
        case PVS:
            for (i = 0; i < lf->v; i++) {
                for (j = 0; j < lf->u; j++) {
                    for (k = 0; k < lf->t; k++) {
                        for (l = 0; l < lf->s; l++) {
                            first = *stream_ptr++;

                            if (depth > 8) {
                                second = *stream_ptr++;

                                lf->val[i][j][k][l] = (img_t) ((first << shift_first) + (second << shift_second));
                                lf->val[i][j][k][l] = (img_t) (lf->val[i][j][k][l] & mask[depth - 8]);
                            }
                            else {
                                lf->val[i][j][k][l] = (img_t) first;
                            }
                        }
                    }
                }
            }

            break;

        case SAI:
        default:
            for (i = 0; i < lf->v; i++) {
                for (k = 0; k < lf->t; k++) {
                    for (j = 0; j < lf->u; j++) {
                        for (l = 0; l < lf->s; l++) {
                            first = *stream_ptr++;

                            if (depth > 8) {
                                second = *stream_ptr++;

                                lf->val[i][j][k][l] = (img_t) ((first << shift_first) + (second << shift_second));
                                lf->val[i][j][k][l] = (img_t) (lf->val[i][j][k][l] & mask[depth - 8]);
                            }
                            else {
                                lf->val[i][j][k][l] = (img_t) first;
                            }
                        }
                    }
                }
            }

            break;
    }

    free(stream);

    fclose(fp);

    return (lf);
}

// TODO: fix functions headers
int ***init_intra_ref_offset(int ts[2], int prd_order) {
    int ***roff, *ptr;
    int t, s, ds, dt, k;

    int min_ds, max_ds, min_dt, max_dt;

    min_ds = max_ds = min_dt = max_dt = 0;

    roff = (int ***) alloc_2d_array(ts[HEIGHT], ts[WIDTH], sizeof(int *));

    //Values to check for special cases
    for (k = 0; k < prd_order; k++) {
        dt = dyx[k].y;
        ds = dyx[k].x;

        if (dt < min_dt) min_dt = dt;
        if (ds < min_ds) min_ds = ds;
        if (dt > max_dt) max_dt = dt;
        if (ds > max_ds) max_ds = ds;
    }

    // Cycle that runs for all the pixels
    for (t = 0; t < ts[HEIGHT]; t++) {
        for (s = 0; s < ts[WIDTH]; s++) {
            ptr = (int *) alloc_mem((prd_order) * sizeof(int));

            //Conditions to check which references are available for each pixel
            if (t == 0) {
                if (s == 0) {
                    roff[t][s] = ptr;
                    ds = 0;
                    dt = ts[HEIGHT];

                    for (k = 0; k < prd_order; k++) {
                        // Points to a line filled with 128 (if max_val = 256)
                        *ptr++ = dt * ts[WIDTH] + ds;
                    }
                }
                else if (s + min_ds <= 0 || s + max_ds >= ts[WIDTH]) {
                    roff[t][s] = ptr;
                    dt = 0;

                    for (k = 0; k < prd_order; k++) {
                        ds = dyx[k].x;

                        if (s + ds < 0) ds = -s;
                        else if (ds >= 0) ds = -1;

                        *ptr++ = dt * ts[WIDTH] + ds;
                    }
                }
                else {
                    roff[t][s] = roff[t][s - 1];
                    free(ptr);
                }
            }
            else if (t + min_dt <= 0) {
                if (s == 0) {
                    roff[t][s] = ptr;

                    for (k = 0; k < prd_order; k++) {
                        dt = dyx[k].y;

                        if (t + dt < 0) dt = -t;
                        else if (dt >= 0) dt = -1;

                        ds = dyx[k].x;

                        if (s + ds < 0) ds = -s;

                        *ptr++ = dt * ts[WIDTH] + ds;
                    }
                }
                else if (s + min_ds <= 0 || s + max_ds >= ts[WIDTH]) {
                    roff[t][s] = ptr;

                    for (k = 0; k < prd_order; k++) {
                        dt = dyx[k].y;

                        if (t + dt < 0) dt = -t;

                        ds = dyx[k].x;

                        if (s + ds < 0) ds = -s;
                        else if (s + ds >= ts[WIDTH]) {
                            ds = ts[WIDTH] - s - 1;
                        }

                        *ptr++ = dt * ts[WIDTH] + ds;
                    }
                }
                else {
                    roff[t][s] = roff[t][s - 1];
                    free(ptr);
                }
            }
            else {
                roff[t][s] = roff[t - 1][s];
                free(ptr);
            }
        }
    }

    return (roff);
}

// TODO: fix functions headers
int ***init_sai_ref_offset(int ts[2], int prd_order, int distance, int ***disparity_vectors, int disp_blk_size) {
    int ***roff, *ptr;
    int dvector[2], disparity = 0;
    int t, s, ds, dt, k, base;

    int min_ds, max_ds, min_dt, max_dt;
    min_ds = max_ds = min_dt = max_dt = 0;

    roff = (int ***) alloc_2d_array(ts[HEIGHT], ts[WIDTH], sizeof(int *));

    //Values to check for special cases
    for (k = 0; k < prd_order - 1; k++) {
        dt = idyx[k].y;
        ds = idyx[k].x;

        if (dt < min_dt) min_dt = dt;
        if (dt > max_dt) max_dt = dt;
        if (ds < min_ds) min_ds = ds;
        if (ds > max_ds) max_ds = ds;
    }

    // TODO: make this better if it works
    // Cycle that runs for all the pixels
    for (int tmp_t = 0; tmp_t < ts[HEIGHT]; tmp_t++) {
        for (int tmp_s = 0; tmp_s < ts[WIDTH]; tmp_s++) {
            ptr = (int *) alloc_mem((prd_order) * sizeof(int));

            if (disparity_vectors != NULL) {
                dvector[0] = disparity_vectors[tmp_t / disp_blk_size][tmp_s / disp_blk_size][0];
                dvector[1] = disparity_vectors[tmp_t / disp_blk_size][tmp_s / disp_blk_size][1];

                t = tmp_t + dvector[0];
                s = tmp_s + dvector[1];
                
                disparity = dvector[1] + (dvector[0] * ts[WIDTH]);
            }
            else {
                t = tmp_t;
                s = tmp_s;
            }

            // Base for reference template placement
            base = ts[WIDTH] * (ts[HEIGHT] + 1) * distance + disparity;

            if (t == 0 || t + min_dt <= 0) {
                if (s == 0) {
                    roff[tmp_t][tmp_s] = ptr;
                    *ptr++ = base;

                    for (k = 0; k < prd_order - 1; k++) {
                        dt = idyx[k].y;
                        ds = idyx[k].x;

                        if (t + dt < 0 || s + ds < 0) {
                            *ptr++ = base;
                        }
                        else {
                            *ptr++ = dt * ts[WIDTH] + ds + base;
                        }
                    }
                }
                else {
                    roff[tmp_t][tmp_s] = ptr;
                    *ptr++ = base;

                    for (k = 0; k < prd_order - 1; k++) {
                        dt = idyx[k].y;
                        ds = idyx[k].x;

                        if (t + dt < 0 || s + ds < 0 || s + ds >= ts[WIDTH]) {
                            *ptr++ = base;
                        }
                        else {
                            *ptr++ = dt * ts[WIDTH] + ds + base;
                        }
                    }
                }
            }
            else {
                roff[tmp_t][tmp_s] = ptr;
                *ptr++ = base;

                for (k = 0; k < prd_order - 1; k++) {
                    dt = idyx[k].y;
                    ds = idyx[k].x;

                    if (t + dt >= ts[HEIGHT] || s + ds < 0 || s + ds >= ts[WIDTH]) {
                        *ptr++ = base;
                    }
                    else {
                        *ptr++ = dt * ts[WIDTH] + ds + base;
                    }
                }
            }
        }
    }

    return (roff);
}

/*--------------------------- init_ref_offset -----------------------*
 |  Function init_ref_offset
 |
 |  Purpose:  Calculates the offset of a reference to a given pixel
 |
 |  Parameters:
 |		img				--> Image structure (IN)
 |		prd_order		--> Number of intra reference pixels (IN)
 |		mi_prd_order	--> Number of neighbour micro images reference pixels (IN)
 |		mi_size     	--> Micro image size (IN)
 |
 |  Returns:  int*** --> Array with the references offset
 *-------------------------------------------------------------------*/
int ***init_ref_offset(int ts[2], int type, int prd_order, int ***disparity, int disp_blk_size) {
    if (type < 0 || type > SAI_REFERENCES) {
        fprintf(stderr, "Wrong usage of init_ref_offset function!\n");
        exit(-10);
    }

    int ***roff;
    
    if (type == 0) {
        roff = init_intra_ref_offset(ts, prd_order);
    }
    else {
        roff = init_sai_ref_offset(ts, prd_order, type, disparity, disp_blk_size);
    }

    return (roff);
}

void free_ref_offset(const int ts[2], int type, int prd_order, int ***roff) {
    if (type < 0 || type > SAI_REFERENCES) {
        fprintf(stderr, "Wrong usage of free_ref_offset function!\n");
        exit(-11);
    }

    int dt, ds, t, s, k;
    int min_ds, max_ds, min_dt, max_dt;
    min_ds = max_ds = min_dt = max_dt = 0;
    
    // Intra
    if (type == 0) {
        // Check if mi_size should be also used
        for (k = 0; k < prd_order; k++) {
            dt = dyx[k].y;
            ds = dyx[k].x;

            if (dt < min_dt) min_dt = dt;
            if (dt > max_dt) max_dt = dt;
            if (ds < min_ds) min_ds = ds;
            if (ds > max_ds) max_ds = ds;
        }

        // Cycle that runs for all the pixels
        for (t = 0; t < ts[HEIGHT]; t++) {
            for (s = 0; s < ts[WIDTH]; s++) {
                //Conditions to check which references are available for each pixel
                if ((t == 0 && (s == 0 || s + min_ds <= 0 || s + max_ds >= ts[WIDTH])) ||
                    (t + min_dt <= 0 && (s == 0 || s + min_ds <= 0 || s + max_ds >= ts[WIDTH]))) {
                        free(roff[t][s]);
                }
            }
        }

        safefree((void **) &(roff));
    }
    else {
        // Check if mi_size should be also used
        for (k = 0; k < prd_order - 1; k++) {
            dt = idyx[k].y;
            ds = idyx[k].x;

            if (dt < min_dt) min_dt = dt;
            if (dt > max_dt) max_dt = dt;
            if (ds < min_ds) min_ds = ds;
            if (ds > max_ds) max_ds = ds;
        }

        //Cycle that runs for all the pixels
        for (t = 0; t < ts[HEIGHT]; t++) {
            for (s = 0; s < ts[WIDTH]; s++) {
                free(roff[t][s]);
            }
        }

        safefree((void **) &(roff));
    }
}

int *init_ctx_weight(int type, int prd_order, int delta) {
    if (type < 0 || type > SAI_REFERENCES) {
        fprintf(stderr, "Wrong usage of free_ref_offset function!\n");
        exit(-11);
    }

    int *ctx_weight, k;
    double dt, ds;

    ctx_weight = (int *) alloc_mem((prd_order) * sizeof(int));

    if (type == 0) {
        for (k = 0; k < prd_order; k++) {
            dt = dyx[k].y;
            ds = dyx[k].x;

            ctx_weight[k] = (int) (64.0 / sqrt(dt * dt + ds * ds) + 0.5);
        }
    }
    else {
        ctx_weight[0] = (int) (64.0 / sqrt(delta * delta) + 0.5);

        for (k = 0; k < prd_order - 1; k++) {
            dt = idyx[k].y;
            ds = idyx[k].x;

            ctx_weight[k + 1] = (int) (64.0 / sqrt(delta * delta + dt * dt + ds * ds) + 0.5);
        }

    }

    return (ctx_weight);
}

int e2E(int e, int prd, int flag, int maxval) {
    int E, th;

    E = (e > 0) ? e : -e;
    th = (prd < ((maxval + 1) >> 1)) ? prd : maxval - prd;

    if (E > th) {
        E += th;
    }
    else if (flag) {
        E = (e < 0) ? (E << 1) - 1 : (E << 1);
    }
    else {
        E = (e > 0) ? (E << 1) - 1 : (E << 1);
    }

    return (E);
}

int E2e(int E, int prd, int flag, int maxval) {
    int e, th;

    th = (prd < ((maxval + 1) >> 1)) ? prd : maxval - prd;

    if (E > (th << 1)) {
        e = (prd < ((maxval + 1) >> 1)) ? E - th : th - E;
    }
    else if (flag) {
        e = (E & 1) ? -((E >> 1) + 1) : (E >> 1);
    }
    else {
        e = (E & 1) ? (E >> 1) + 1 : -(E >> 1);
    }

    return (e);
}

void mtf_classlabel(char **class, int *mtfbuf, int t, int s, int bsize, int width, int num_class) {
    int i, j, k, ref[3];

    if (t == 0) {
        if (s == 0) {
            ref[0] = ref[1] = ref[2] = 0;
        }
        else {
            ref[0] = ref[1] = ref[2] = (int) class[t][s - 1];
        }
    }
    else {
        ref[0] = (int) class[t - 1][s];
        ref[1] = (s == 0) ? class[t - 1][s] : class[t][s - 1];
        ref[2] = (s + bsize >= width) ? class[t - 1][s] : class[t - 1][s + bsize];
        if (ref[1] == ref[2]) {
            ref[2] = ref[0];
            ref[0] = ref[1];
        }
    }

    /* move to front */
    for (k = 2; k >= 0; k--) {
        if ((j = mtfbuf[ref[k]]) == 0) continue;

        for (i = 0; i < num_class; i++) {
            if (mtfbuf[i] < j) {
                mtfbuf[i]++;
            }
        }

        mtfbuf[ref[k]] = 0;
    }
}

// Returns time
double cpu_time(void) {
#include <time.h>
#ifndef HAVE_CLOCK
    #  include <sys/times.h>
	struct tms t;
#endif
#ifndef CLK_TCK
#  define CLK_TCK 60
#endif
    static clock_t prev = 0;
    clock_t cur, dif;

#ifdef HAVE_CLOCK
    cur = clock();
#else
    times(&t);
	cur = t.tms_utime + t.tms_stime;
#endif
    if (cur > prev) {
        dif = cur - prev;
    }
    else {
        dif = (unsigned) cur - prev;
    }
    prev = cur;

#ifdef HAVE_CLOCK
    return ((double) dif / CLOCKS_PER_SEC);
#else
    return ((double)dif / CLK_TCK);
#endif
}

char* cat_str(char* str1, char* str2, int type) {
    char *new_str = (char *) alloc_mem(sizeof(char) * (strlen(str1) + strlen(str2) + 1));

    if(new_str != NULL){
        new_str[0] = '\0';   // ensures the memory is an empty string

        strcat(new_str, str1);
        strcat(new_str, str2);
    }
    else {
        printf("malloc failed!\n");
        exit(-12);
    }

    if (type == 1 ) safefree((void **) &str1);
    safefree((void **) &str2);

    return new_str;
}

// buffer must have length >= sizeof(int) + 1
// Write to the buffer backwards so that the binary representation
// is in the correct order i.e.  the LSB is on the far right
// instead of the far left of the printed string
char *int2bin(int a, int buf_size) {
    char *buffer = (char *)alloc_mem(sizeof(char) * (buf_size + 1));
    buffer[buf_size] = '\0';

    for (int i = buf_size - 1; i >= 0; i--) {
        buffer[i] = (char) ((a & 1) + '0');

        a >>= 1;
    }

    return buffer;
}

// Convert the frame number to coordinates in VxU
void frame2coordinates(int *coordinates, int frame, const int U) {
    coordinates[0] = frame / U;
    coordinates[1] = frame % U;
}

// Compare SAI distance for the sort function
int compare_distance(const void *a, const void *b) {
    SAIDISTANCE *saidist_a = *((SAIDISTANCE**) a);
    SAIDISTANCE *saidist_b = *((SAIDISTANCE**) b);

    double diff = saidist_a->distance - saidist_b->distance;
    int signal = diff < 0 ? -1 : diff > 0 ? 1 : 0;

    return (signal);
}

// Function for calculating median as seen in: https://en.wikiversity.org/wiki/C_Source_Code/Find_the_median_and_mean
double median(int vector[], int n) {
    int temp;
    int i, j;

    // the following two loops sort the array vector in ascending order
    for(i = 0; i < n - 1; i++) {
        for(j = i + 1; j < n; j++) {
            if(vector[j] < vector[i]) {
                // swap elements
                temp = vector[i];
                vector[i] = vector[j];
                vector[j] = temp;
            }
        }
    }

    if(n % 2 == 0) {
        // if there is an even number of elements, return mean of the two elements in the middle
        return((vector[n / 2] + vector[n / 2 - 1]) / 2.0);
    } else {
        // else return the element in the middle
        return vector[n / 2];
    }
}

int get_random_access_region(int no_ra_regions, const int vu[2], int v, int u, int vu_curr[2]) {
    int ra_regions2[13][13] = {
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
            { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
            { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
            { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
            { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
            { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
            { 1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1,  1},
    };

    int ra_regions4[13][13] = {
            { 0,  0,  0,  0,  0,  0, -1,  1,  1,  1,  1,  1,  1},
            { 0,  0,  0,  0,  0,  0, -1,  1,  1,  1,  1,  1,  1},
            { 0,  0,  0,  0,  0,  0, -1,  1,  1,  1,  1,  1,  1},
            { 0,  0,  0,  0,  0,  0, -1,  1,  1,  1,  1,  1,  1},
            { 0,  0,  0,  0,  0,  0, -1,  1,  1,  1,  1,  1,  1},
            { 0,  0,  0,  0,  0,  0, -1,  1,  1,  1,  1,  1,  1},
            {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1},
            { 2,  2,  2,  2,  2,  2, -1,  3,  3,  3,  3,  3,  3},
            { 2,  2,  2,  2,  2,  2, -1,  3,  3,  3,  3,  3,  3},
            { 2,  2,  2,  2,  2,  2, -1,  3,  3,  3,  3,  3,  3},
            { 2,  2,  2,  2,  2,  2, -1,  3,  3,  3,  3,  3,  3},
            { 2,  2,  2,  2,  2,  2, -1,  3,  3,  3,  3,  3,  3},
            { 2,  2,  2,  2,  2,  2, -1,  3,  3,  3,  3,  3,  3},
    };

    int ra_regions5[13][13] = {
            { 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0},
            { 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  3},
            { 1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  3,  3},
            { 1,  1,  1,  0,  0,  0,  0,  0,  0,  0,  3,  3,  3},
            { 1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3},
            { 1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3},
            { 1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3},
            { 1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3},
            { 1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3,  3},
            { 1,  1,  1,  4,  4,  4,  4,  4,  4,  4,  3,  3,  3},
            { 1,  1,  4,  4,  4,  4,  4,  4,  4,  4,  3,  3,  3},
            { 1,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  3},
            { 4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4,  4},
    };

    int ra_regions9[13][13] = {
            { 0,  0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  2,  2},
            { 0,  0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  2,  2},
            { 0,  0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  2,  2},
            { 0,  0,  0,  0,  1,  1,  1,  1,  1,  2,  2,  2,  2},
            { 3,  3,  3,  3,  4,  4,  4,  4,  4,  5,  5,  5,  5},
            { 3,  3,  3,  3,  4,  4,  4,  4,  4,  5,  5,  5,  5},
            { 3,  3,  3,  3,  4,  4,  4,  4,  4,  5,  5,  5,  5},
            { 3,  3,  3,  3,  4,  4,  4,  4,  4,  5,  5,  5,  5},
            { 3,  3,  3,  3,  4,  4,  4,  4,  4,  5,  5,  5,  5},
            { 6,  6,  6,  6,  7,  7,  7,  7,  7,  8,  8,  8,  8},
            { 6,  6,  6,  6,  7,  7,  7,  7,  7,  8,  8,  8,  8},
            { 6,  6,  6,  6,  7,  7,  7,  7,  7,  8,  8,  8,  8},
            { 6,  6,  6,  6,  7,  7,  7,  7,  7,  8,  8,  8,  8},
    };

    int vu_stride[2];
    int ra_region, curr_ra_region;

    vu_stride[HEIGHT] = (vu[HEIGHT] == 13) ? 0 : 2;
    vu_stride[WIDTH] = (vu[WIDTH] == 13) ? 0 : 2;

    switch (no_ra_regions) {
        case 2:
            ra_region = ra_regions2[v + vu_stride[HEIGHT]][u + vu_stride[WIDTH]];
            break;

        case 4:
            ra_region = ra_regions4[v + vu_stride[HEIGHT]][u + vu_stride[WIDTH]];
            break;

        case 5:
            ra_region = ra_regions5[v + vu_stride[HEIGHT]][u + vu_stride[WIDTH]];
            break;

        case 9:
            ra_region = ra_regions9[v + vu_stride[HEIGHT]][u + vu_stride[WIDTH]];
            break;

        default:
            fprintf(stderr, "Something went wrong. Wrong number of Random Access Regions: %d.\n", no_ra_regions);
            exit(-9);
    }

    if (vu_curr != NULL) {
        curr_ra_region = get_random_access_region(no_ra_regions, vu, vu_curr[HEIGHT], vu_curr[WIDTH], NULL);

        if (ra_region != -1 && curr_ra_region == -1) {
            if (no_ra_regions == 2) {
                ra_region = -1;
            }
            else {
                if ((vu_curr[HEIGHT] == vu[HEIGHT] / 2 &&
                     ((vu_curr[WIDTH] <= vu[WIDTH] / 2 && u < vu[WIDTH] / 2) ||
                      (vu_curr[WIDTH] >= vu[WIDTH] / 2 && u > vu[WIDTH] / 2))) ||
                    ((vu_curr[WIDTH] == vu[WIDTH] / 2 &&
                      ((vu_curr[HEIGHT] <= vu[HEIGHT] / 2 && u < vu[HEIGHT] / 2) ||
                       (vu_curr[HEIGHT] >= vu[HEIGHT] / 2 && u > vu[HEIGHT] / 2)))) ||
                    ((vu_curr[HEIGHT] == vu[HEIGHT] / 2 && vu_curr[WIDTH] == vu[WIDTH] / 2))) {
                    ra_region = -1;
                }
            }
        }
        else if (ra_region == -1) {
            if (curr_ra_region == -1) {
                if (no_ra_regions == 4 && ((vu_curr[HEIGHT] > vu[HEIGHT] / 2 && v < vu[HEIGHT] / 2) ||
                                           (vu_curr[HEIGHT] < vu[HEIGHT] / 2 && v > vu[HEIGHT] / 2))) {
                    ra_region = 10;
                }
            }
            else {
                if (no_ra_regions == 2) {
                    ra_region = curr_ra_region;
                }
                else {
                    if (curr_ra_region == 0 && ((v == vu[HEIGHT] / 2 && vu_curr[WIDTH] <= vu[WIDTH] / 2) ||
                                                (v == vu[WIDTH] / 2 && vu_curr[HEIGHT] <= vu[HEIGHT] / 2))) {
                        ra_region = 0;
                    }
                    else if (curr_ra_region == 1 && ((v == vu[HEIGHT] / 2 && vu_curr[WIDTH] >= vu[WIDTH] / 2) ||
                                                     (v == vu[WIDTH] / 2 && vu_curr[HEIGHT] <= vu[HEIGHT] / 2))) {
                        ra_region = 1;
                    }
                    else if (curr_ra_region == 2 && ((v == vu[HEIGHT] / 2 && vu_curr[WIDTH] <= vu[WIDTH] / 2) ||
                                                     (v == vu[WIDTH] / 2 && vu_curr[HEIGHT] >= vu[HEIGHT] / 2))) {
                        ra_region = 2;
                    }
                    else if (curr_ra_region == 3 && ((v == vu[HEIGHT] / 2 && vu_curr[WIDTH] >= vu[WIDTH] / 2) ||
                                                     (v == vu[WIDTH] / 2 && vu_curr[HEIGHT] >= vu[HEIGHT] / 2))) {
                        ra_region = 3;
                    }
                }
            }
        }
    }

    return ra_region;
}