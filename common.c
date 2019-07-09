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

//Reference pixel position separated by their distances
const POINT tridyx[] = {
        /* 0 --> 1 */
        //{ 0, 0},
        /* 1 --> 4 */
        { 0,-1}, {-1, 0}, { 0, 1},
        /* 2 --> 9 */
        { 0,-2}, {-1,-1}, {-2, 0}, {-1, 1}, { 0, 2},
        /* 3 --> 16 */
        { 0,-3}, {-1,-2}, {-2,-1}, {-3, 0}, {-2, 1}, {-1, 2}, { 0, 3},
        /* 4 --> 25 */
        { 0,-4}, {-1,-3}, {-2,-2}, {-3,-1}, {-4, 0}, {-3, 1}, {-2, 2}, {-1, 3}, { 0, 4},
        /* 5 --> 36 */
        { 0,-5}, {-1,-4}, {-2,-3}, {-3,-2}, {-4,-1}, {-5, 0}, {-4, 1}, {-3, 2}, {-2, 3}, {-1, 4}, { 0, 5},
        /* 6 --> 49 */
        { 0,-6}, {-1,-5}, {-2,-4}, {-3,-3}, {-4,-2}, {-5,-1}, {-6, 0}, {-5, 1}, {-4, 2}, {-3, 3}, {-2, 4}, {-1, 5}, { 0, 6},
};

int bref1[2][5] = {{2, -2, -1, 0, -1},
                   {-1, -1, 1, 0, 2},
};

int bref2[3][5] = {{3, -3, -1, 0, -1},
                   {-2, -1, 2, 0, 3},
                   {1, -1, 1, 1, 3}
};

int bref3[4][5] = {{4, -4, -1, 0, -1},
                   {-2, -2, 2, 0, 4},
                   {-1, -1, 1, 0, 2},
                   {2, -1, 1, 2, 4}
};

int bref4[5][5] = {{5, -5, -1, 0, -1},
                   {-3, -2, 3, 0, 5},
                   {-1, -1, 1, 0, 2},
                   {2, -1, 2, 2, 5},
                   {1, -1, 1, 3, 5}
};

int bref5[6][5] = {{6, -6, -1, 0, -1},
                   {-3, -3, 3, 0, 6},
                   {-2, -1, 2, 0, 3},
                   {1, -1, 1, 1, 3},
                   {2, -1, 2, 3, 6},
                   {1, -1, 1, 4, 6}
};

int bref6[7][5] = {{7, -7, -1, 0, -1},
                   {-4, -3, 4, 0, 7},
                   {-2, -1, 2, 0, 3},
                   {1, -1, 1, 1, 3},
                   {3, -2, 2, 3, 7},
                   {-1, -1, 1, 3, 5},
                   {2, -1, 1, 5, 7}
};

int bref7[8][5] = {{8, -8, -1, 0, -1},
                   {-4, -4, 4, 0, 8},
                   {-2, -2, 2, 0, 4},
                   {-1, -1, 1, 0, 2},
                   {2, -1, 1, 2, 4},
                   {3, -2, 2, 4, 8},
                   {-1, -1, 1, 4, 8},
                   {2, -1, 1, 6, 8}
};

int bref8[9][5] = {{9, -9, -1, 0, -1},
                   {-5, -4, 5, 0, 9},
                   {-2, -2, 2, 0, 4},
                   {-1, -1, 1, 0, 2},
                   {2, -1, 1, 2, 4},
                   {3, -2, 3, 4, 9},
                   {-1, -1, 1, 4, 6},
                   {2, -1, 2, 6, 9},
                   {1, -1, 1, 7, 9}
};

int bref9[10][5] = {{10, -10, -1, 0, -1},
                    {-5, -5, 5, 0, 10},
                    {-3, -2, 3, 0, 5},
                    {-1, -1, 1, 0, 2},
                    {2, -1, 2, 2, 5},
                    {1, -1, 1, 3, 5},
                    {3, -2, 3, 5, 10},
                    {-1, -1, 1, 5, 7},
                    {2, -1, 2, 7, 10},
                    {1, -1, 1, 8, 10}
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
void ***alloc_3d_array(size_t height, size_t width, size_t frames, size_t size) {
    void ***mat;
    int i, j;

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
    int i, j, k;

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

    /*for (i = 0; i < v; i++) {
        mat[i] = &vv[i * u];
    }

    for (i = 0; i < v; i++) {
        for (j = 0; j < u; j++) {
            vv[i * u + j] = &vu[i * u * t * size + j * t * size];
        }
    }

    for (i = 0; i < v; i++) {
        for (j = 0; j < u; j++) {
            for (k = 0; k < t; k++) {
                vu[i * u + j * t + k * size] = &vut[i * u * t * s * size + j * t * s * size + k * s * size];
            }
        }
    }*/

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

// Function to copy LF4D
LF4D *copy_yuv(LF4D *lf) {
    int i, j, k, l;

    LF4D *new_lf = alloc_lf4d(lf->v, lf->u, lf->t, lf->s, lf->maxval);

    for (i = 0; i < lf->v; i++) {
        for (j = 0; j < lf->u; j++) {
            for (k = 0; k < lf->t; k++) {
                for (l = 0; l < lf->s; l++) {
                    new_lf->val[i][j][k][l] = lf->val[i][j][k][l];
                }
            }
        }
    }

    return (new_lf);
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

//    int dimensions[2];
//    int y
//    int x;
//    int* idx1;
//    int* idx2;
//
//    if( loopyfirst )
//    {
//        idx0 = &y;
//        idx1 = &x;
//        dimensions[0] = height
//        dimensions[1] = width
//    }
//    else
//    {
//        idx0 = &x;
//        idx1 = &y;
//        dimensions[0] = width
//        dimensions[1] = height
//    }
//
//    while( *idx1 < dimensions[0] )
//    {
//        while( *idx2 < dimensions[1] )
//        {
//            img[y][x] = input;
//            *idx2++;
//        }
//        *idx1++;
//    }

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
 |		depth		--> Image dynamic range (bpp) (IN)
 |		endianness	--> YUV endianness, for depth > 8 bpp (IN)
 |		format      --> Light field file format (IN)
 |
 |  Returns:  IMAGE* --> returns a video type structure
 *-------------------------------------------------------------------*/
LF4D *read_yuv(char *filename, int v, int u, int t, int s, int depth, int endianness, int chroma, int format) {
    int i, j, k, l;
    unsigned int shift_first, shift_second, first, second;
    double chroma_pass = (chroma == GRAY ? 1 : (chroma == S444 ? 3 : 1.5));
    LF4D *lf;
    FILE *fp;


    int elements = v * u * t * s * (depth == 8 ? 1 : 2);
    unsigned char *stream_ptr, *stream = (unsigned char *) alloc_mem(elements * sizeof(unsigned char));

    //Open file
    fp = fileopen(filename, "rb");

    // Check if image dimensions are correct (It has to be multiple of BASE_BSIZE)
//    if ((width % BASE_BSIZE) || (height % BASE_BSIZE)) {
//        fprintf(stderr, "Image width and height must be multiples of %d!\n", BASE_BSIZE);
//        exit(1);
//    }

    fread(stream, sizeof(unsigned char), elements, fp);
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

    fclose(fp);

    return (lf);
}

// TODO: fix functions headers
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
int *****init_ref_offset(int vu[2], int ts[2], int type, int prd_order) {
    int *****roff, *ptr;
    int v, u, t, s, ds, dt, k, base = 0;

    int min_ds, max_ds, min_dt, max_dt;
    int mmin_ds, mmax_ds, mmin_dt, mmax_dt;

    min_ds = max_ds = min_dt = max_dt = 0;
    mmin_ds = mmax_ds = mmin_dt = mmax_dt = 0;

    roff = (int *****) alloc_4d_array(vu[HEIGHT], vu[WIDTH], ts[HEIGHT], ts[WIDTH], sizeof(int *));

    switch (type) {
        case INTRA_PRED: // Intra
        case MI_BORDERS_PRED:
            base = 0;
            break;

        case MI_UP_PRED:
            base = -ts[WIDTH] * ts[HEIGHT] * vu[WIDTH];
            break;

        case MI_LEFT_PRED:
            base = -ts[WIDTH] * ts[HEIGHT];
            break;

        case MI_LDIAG_PRED:
            base = -ts[WIDTH] * ts[HEIGHT] * (vu[WIDTH] + 1);
            break;

        case MI_RDIAG_PRED:
            base = -ts[WIDTH] * ts[HEIGHT] * (vu[WIDTH] - 1);
            break;

        default:
            fprintf(stderr, "Wrong usage of init_ref_offset function!\n");
            exit(-10);
    }

    switch (type) {
        case INTRA_PRED: // Intra
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
            for (v = 0; v < vu[HEIGHT]; v++) {
                for (u = 0; u < vu[WIDTH]; u++) {
                    for (t = 0; t < ts[HEIGHT]; t++) {
                        for (s = 0; s < ts[WIDTH]; s++) {
                            ptr = (int *) alloc_mem((prd_order) * sizeof(int));

                            //Conditions to check which references are available for each pixel
                            if (v == 0 && u == 0) {
                                if (t == 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
                                        ds = 0;
                                        dt = ts[HEIGHT];

                                        for (k = 0; k < prd_order; k++) {
                                            // Points to a line filled with 128 (if max_val = 256)
                                            *ptr++ = dt * ts[WIDTH] * (vu[HEIGHT] * vu[WIDTH] - v * vu[WIDTH] - u) + ds;
                                        }
                                    }
                                    else if (s + min_ds <= 0 || s + max_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
                                        dt = 0;

                                        for (k = 0; k < prd_order; k++) {
                                            ds = dyx[k].x;

                                            if (s + ds < 0) ds = -s;
                                            else if (ds >= 0) ds = -1;

                                            *ptr++ = dt * ts[WIDTH] + ds;
                                        }
                                    }
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else if (t + min_dt <= 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;

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
                                        roff[v][u][t][s] = ptr;

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
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u][t - 1][s];
                                    free(ptr);
                                }
                            }
                            else if (u > 0) {
                                if (t == 0 && s == 0) {
                                    roff[v][u][t][s] = ptr;
                                    ds = 0;
                                    dt = ts[HEIGHT];

                                    for (k = 0; k < prd_order; k++) {
                                        // Points to a line filled with 128 (if max_val = 256)
                                        *ptr++ = dt * ts[WIDTH] * (vu[HEIGHT] * vu[WIDTH] - v * vu[WIDTH] - u) + ds;
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u - 1][t][s];
                                    free(ptr);
                                }
                            }
                            else {
                                if (t == 0 && s == 0) {
                                    roff[v][u][t][s] = ptr;
                                    ds = 0;
                                    dt = ts[HEIGHT];

                                    for (k = 0; k < prd_order; k++) {
                                        // Points to a line filled with 128 (if max_val = 256)
                                        *ptr++ = dt * ts[WIDTH] * (vu[HEIGHT] * vu[WIDTH] - v * vu[WIDTH] - u) + ds;
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v - 1][u][t][s];
                                    free(ptr);
                                }
                            }
                        }
                    }
                }
            }

            break;

        case MI_UP_PRED: // Up prediction
            for (k = 0; k < prd_order; k++) {
                dt = dyx[k].y;
                ds = dyx[k].x;

                if (dt < min_dt) min_dt = dt;
                if (ds < min_ds) min_ds = ds;
                if (ds > max_ds) max_ds = ds;
            }

            for (k = 0; k < prd_order - 1; k++) {
                dt = idyx[k].y;
                ds = idyx[k].x;

                if (dt < mmin_dt) mmin_dt = dt;
                if (dt > mmax_dt) mmax_dt = dt;
                if (ds < mmin_ds) mmin_ds = ds;
                if (ds > mmax_ds) mmax_ds = ds;
            }

            // Cycle that runs for all the pixels
            for (v = 0; v < vu[HEIGHT]; v++) {
                for (u = 0; u < vu[WIDTH]; u++) {
                    for (t = 0; t < ts[HEIGHT]; t++) {
                        for (s = 0; s < ts[WIDTH]; s++) {
                            ptr = (int *) alloc_mem((prd_order) * sizeof(int));

                            // LF reference offset
                            if (v == 0 && u == 0) {
                                if (t == 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
                                        ds = 0;
                                        dt = ts[HEIGHT];

                                        for (k = 0; k < prd_order; k++) {
                                            // Points to a line filled with 128 (if max_val = 256)
                                            *ptr++ = dt * ts[WIDTH] * (vu[HEIGHT] * vu[WIDTH] - v * vu[WIDTH] - u) + ds;
                                        }
                                    }
                                    else if (s + min_ds <= 0 || s + max_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
                                        dt = 0;

                                        for (k = 0; k < prd_order; k++) {
                                            ds = dyx[k].x;

                                            if (s + ds < 0) ds = -s;
                                            else if (ds >= 0) ds = -1;

                                            *ptr++ = dt * ts[WIDTH] + ds;
                                        }
                                    }
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else if (t + min_dt <= 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;

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
                                        roff[v][u][t][s] = ptr;

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
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u][t - 1][s];
                                    free(ptr);
                                }
                            }
                            else if (v == 0) {
                                if (t == 0 && s == 0) {
                                    roff[v][u][t][s] = ptr;
                                    ds = 0;
                                    dt = ts[HEIGHT];

                                    for (k = 0; k < prd_order; k++) {
                                        // Points to a line filled with 128 (if max_val = 256)
                                        *ptr++ = dt * ts[WIDTH] * (vu[HEIGHT] * vu[WIDTH] - v * vu[WIDTH] - u) + ds;
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u - 1][t][s];
                                    free(ptr);
                                }
                            }
                            else if (v == 1 && u == 0) {
                                if (t == 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
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
                                    else if (s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
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
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else if (t + mmin_dt <= 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
                                        *ptr++ = base;

                                        for (k = 0; k < prd_order - 1; k++) {
                                            dt = idyx[k].y;
                                            ds = idyx[k].x;

                                            if (t + dt < ts[HEIGHT] || s + ds < 0) {
                                                *ptr++ = base;
                                            }
                                            else {
                                                *ptr++ = dt * ts[WIDTH] + ds + base;
                                            }
                                        }
                                    }
                                    else if (s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
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
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else if (t + mmax_dt >= ts[HEIGHT]) {
                                    if (s == 0 || s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
                                        *ptr++ = base;

                                        for (k = 0; k < prd_order - 1; k++) {
                                            dt = idyx[k].y;
                                            ds = idyx[k].x;

                                            if (t + dt >= ts[HEIGHT] || s + ds < 0 || s + ds >= ts[HEIGHT]) {
                                                *ptr++ = base;
                                            }
                                            else {
                                                *ptr++ = dt * ts[WIDTH] + ds + base;
                                            }
                                        }
                                    }
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u][t - 1][s];
                                    free(ptr);
                                }
                            }
                            else if (v == 1) {
                                roff[v][u][t][s] = roff[v][u - 1][t][s];
                                free(ptr);
                            }
                            else {
                                roff[v][u][t][s] = roff[v - 1][u][t][s];
                                free(ptr);
                            }
                        }
                    }
                }
            }

            break;

        case MI_LEFT_PRED: // Left prediction
            for (k = 0; k < prd_order; k++) {
                dt = dyx[k].y;
                ds = dyx[k].x;

                if (dt < min_dt) min_dt = dt;
                if (ds < min_ds) min_ds = ds;
                if (ds > max_ds) max_ds = ds;
            }

            for (k = 0; k < prd_order - 1; k++) {
                dt = tridyx[k].y;
                ds = tridyx[k].x;

                if (dt < mmin_dt) mmin_dt = dt;
                if (dt > mmax_dt) mmax_dt = dt;
                if (ds < mmin_ds) mmin_ds = ds;
                if (ds > mmax_ds) mmax_ds = ds;
            }

            // Cycle that runs for all the pixels
            for (v = 0; v < vu[HEIGHT]; v++) {
                for (u = 0; u < vu[WIDTH]; u++) {
                    for (t = 0; t < ts[HEIGHT]; t++) {
                        for (s = 0; s < ts[WIDTH]; s++) {
                            ptr = (int *) alloc_mem((prd_order) * sizeof(int));

                            // LF reference offset
                            if (v == 0 && u == 0) {
                                if (t == 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
                                        ds = 0;
                                        dt = ts[HEIGHT];

                                        for (k = 0; k < prd_order; k++) {
                                            // Points to a line filled with 128 (if max_val = 256)
                                            *ptr++ = dt * ts[WIDTH] * (vu[HEIGHT] * vu[WIDTH] - v * vu[WIDTH] - u) + ds;
                                        }
                                    }
                                    else if (s + min_ds <= 0 || s + max_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
                                        dt = 0;

                                        for (k = 0; k < prd_order; k++) {
                                            ds = dyx[k].x;

                                            if (s + ds < 0) ds = -s;
                                            else if (ds >= 0) ds = -1;

                                            *ptr++ = dt * ts[WIDTH] + ds;
                                        }
                                    }
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else if (t + min_dt <= 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;

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
                                        roff[v][u][t][s] = ptr;

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
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u][t - 1][s];
                                    free(ptr);
                                }
                            }
                            else if (u == 0) {
                                if (t == 0 && s == 0) {
                                    roff[v][u][t][s] = ptr;
                                    ds = 0;
                                    dt = ts[HEIGHT];

                                    for (k = 0; k < prd_order; k++) {
                                        // Points to a line filled with 128 (if max_val = 256)
                                        *ptr++ = dt * ts[WIDTH] * (vu[HEIGHT] * vu[WIDTH] - v * vu[WIDTH] - u) + ds;
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v - 1][u][t][s];
                                    free(ptr);
                                }
                            }
                            else if (v == 0 && u == 1) {
                                if (t == 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
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
                                    else if (s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
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
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else if (t + mmin_dt <= 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
                                        *ptr++ = base;

                                        for (k = 0; k < prd_order - 1; k++) {
                                            dt = tridyx[k].y;
                                            ds = tridyx[k].x;

                                            if (t + dt < 0 || s + ds < 0) {
                                                *ptr++ = base;
                                            }
                                            else {
                                                *ptr++ = dt * ts[WIDTH] + ds + base;
                                            }
                                        }
                                    }
                                    else if (s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
                                        *ptr++ = base;

                                        for (k = 0; k < prd_order - 1; k++) {
                                            dt = tridyx[k].y;
                                            ds = tridyx[k].x;

                                            if (t + dt < 0 || s + ds < 0 || s + ds >= ts[WIDTH]) {
                                                *ptr++ = base;
                                            }
                                            else {
                                                *ptr++ = dt * ts[WIDTH] + ds + base;
                                            }
                                        }
                                    }
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u][t - 1][s];
                                    free(ptr);
                                }
                            }
                            else if (u == 1){
                                roff[v][u][t][s] = roff[v - 1][u][t][s];
                                free(ptr);
                            }
                            else {
                                roff[v][u][t][s] = roff[v][u - 1][t][s];
                                free(ptr);
                            }
                        }
                    }
                }
            }

            break;

        case MI_LDIAG_PRED: // Left diagonal prediction
            for (k = 0; k < prd_order; k++) {
                dt = dyx[k].y;
                ds = dyx[k].x;

                if (dt < min_dt) min_dt = dt;
                if (ds < min_ds) min_ds = ds;
                if (ds > max_ds) max_ds = ds;
            }

            for (k = 0; k < prd_order - 1; k++) {
                dt = idyx[k].y;
                ds = idyx[k].x;

                if (dt < mmin_dt) mmin_dt = dt;
                if (dt > mmax_dt) mmax_dt = dt;
                if (ds < mmin_ds) mmin_ds = ds;
                if (ds > mmax_ds) mmax_ds = ds;
            }

            // Cycle that runs for all the pixels
            for (v = 0; v < vu[HEIGHT]; v++) {
                for (u = 0; u < vu[WIDTH]; u++) {
                    for (t = 0; t < ts[HEIGHT]; t++) {
                        for (s = 0; s < ts[WIDTH]; s++) {
                            ptr = (int *) alloc_mem((prd_order) * sizeof(int));

                            // LF reference offset
                            if (v == 0 || u == 0) {
                                if (t == 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
                                        ds = 0;
                                        dt = ts[HEIGHT];

                                        for (k = 0; k < prd_order; k++) {
                                            // Points to a line filled with 128 (if max_val = 256)
                                            *ptr++ = dt * ts[WIDTH] * (vu[HEIGHT] * vu[WIDTH] - v * vu[WIDTH] - u) + ds;
                                        }
                                    }
                                    else if (s + min_ds <= 0 || s + max_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
                                        dt = 0;

                                        for (k = 0; k < prd_order; k++) {
                                            ds = dyx[k].x;

                                            if (s + ds < 0) ds = -s;
                                            else if (ds >= 0) ds = -1;

                                            *ptr++ = dt * ts[WIDTH] + ds;
                                        }
                                    }
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else if (t + min_dt <= 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;

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
                                        roff[v][u][t][s] = ptr;

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
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u][t - 1][s];
                                    free(ptr);
                                }
                            }
                            else if (u == 0) {
                                if (t == 0 && s == 0) {
                                    roff[v][u][t][s] = ptr;
                                    ds = 0;
                                    dt = ts[HEIGHT];

                                    for (k = 0; k < prd_order; k++) {
                                        // Points to a line filled with 128 (if max_val = 256)
                                        *ptr++ = dt * ts[WIDTH] * (vu[HEIGHT] * vu[WIDTH] - v * vu[WIDTH] - u) + ds;
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v - 1][u][t][s];
                                    free(ptr);
                                }
                            }
                            else if (v == 0) {
                                if (t == 0 && s == 0) {
                                    roff[v][u][t][s] = ptr;
                                    ds = 0;
                                    dt = ts[HEIGHT];

                                    for (k = 0; k < prd_order; k++) {
                                        // Points to a line filled with 128 (if max_val = 256)
                                        *ptr++ = dt * ts[WIDTH] * (vu[HEIGHT] * vu[WIDTH] - v * vu[WIDTH] - u) + ds;
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u - 1][t][s];
                                    free(ptr);
                                }
                            }
                            else if (v == 1 && u == 1) {
                                if (t == 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
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
                                    else if (s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
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
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else if (t + mmin_dt <= 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
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
                                    else if (s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
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
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else if (t + mmax_dt >= ts[HEIGHT]) {
                                    if (s == 0 || s + mmin_ds <= 0 || s + max_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
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
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u][t - 1][s];
                                    free(ptr);
                                }
                            }
                            else if (v == 1) {
                                roff[v][u][t][s] = roff[v][u - 1][t][s];
                                free(ptr);
                            }
                            else {
                                roff[v][u][t][s] = roff[v - 1][u][t][s];
                                free(ptr);
                            }
                        }
                    }
                }
            }

            break;

        case MI_RDIAG_PRED: // Right diagonal prediction
            for (k = 0; k < prd_order; k++) {
                dt = dyx[k].y;
                ds = dyx[k].x;

                if (dt < min_dt) min_dt = dt;
                if (ds < min_ds) min_ds = ds;
                if (ds > max_ds) max_ds = ds;
            }

            for (k = 0; k < prd_order - 1; k++) {
                dt = idyx[k].y;
                ds = idyx[k].x;

                if (dt < mmin_dt) mmin_dt = dt;
                if (dt > mmax_dt) mmax_dt = dt;
                if (ds < mmin_ds) mmin_ds = ds;
                if (ds > mmax_ds) mmax_ds = ds;
            }

            // Cycle that runs for all the pixels
            for (v = 0; v < vu[HEIGHT]; v++) {
                for (u = 0; u < vu[WIDTH]; u++) {
                    for (t = 0; t < ts[HEIGHT]; t++) {
                        for (s = 0; s < ts[WIDTH]; s++) {
                            ptr = (int *) alloc_mem((prd_order) * sizeof(int));

                            // LF reference offset
                            if (v == 0 || u == vu[WIDTH] - 1) {
                                if (t == 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
                                        ds = 0;
                                        dt = ts[HEIGHT];

                                        for (k = 0; k < prd_order; k++) {
                                            // Points to a line filled with 128 (if max_val = 256)
                                            *ptr++ = dt * ts[WIDTH] * (vu[HEIGHT] * vu[WIDTH] - v * vu[WIDTH] - u) + ds;
                                        }
                                    }
                                    else if (s + min_ds <= 0 || s + max_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
                                        dt = 0;

                                        for (k = 0; k < prd_order; k++) {
                                            ds = dyx[k].x;

                                            if (s + ds < 0) ds = -s;
                                            else if (ds >= 0) ds = -1;

                                            *ptr++ = dt * ts[WIDTH] + ds;
                                        }
                                    }
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else if (t + min_dt <= 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;

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
                                        roff[v][u][t][s] = ptr;

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
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u][t - 1][s];
                                    free(ptr);
                                }
                            }
                            else if (v == 0) {
                                if (t == 0 && s == 0) {
                                    roff[v][u][t][s] = ptr;
                                    ds = 0;
                                    dt = ts[HEIGHT];

                                    for (k = 0; k < prd_order; k++) {
                                        // Points to a line filled with 128 (if max_val = 256)
                                        *ptr++ = dt * ts[WIDTH] * (vu[HEIGHT] * vu[WIDTH] - v * vu[WIDTH] - u) + ds;
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u - 1][t][s];
                                    free(ptr);
                                }
                            }
                            else if (u == vu[WIDTH] - 1) {
                                roff[v][u][t][s] = roff[v - 1][u][t][s];
                                free(ptr);
                            }
                            else if (v == 1 && u == 0) {
                                if (t == 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
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
                                    else if (s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
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
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else if (t + mmin_dt <= 0) {
                                    if (s == 0) {
                                        roff[v][u][t][s] = ptr;
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
                                    else if (s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
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
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else if (t + mmax_dt >= ts[HEIGHT]) {
                                    if (s == 0 || s + mmin_ds <= 0 || s + max_ds >= ts[WIDTH]) {
                                        roff[v][u][t][s] = ptr;
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
                                    else {
                                        roff[v][u][t][s] = roff[v][u][t][s - 1];
                                        free(ptr);
                                    }
                                }
                                else {
                                    roff[v][u][t][s] = roff[v][u][t - 1][s];
                                    free(ptr);
                                }
                            }
                            else if (v == 1) {
                                roff[v][u][t][s] = roff[v][u - 1][t][s];
                                free(ptr);
                            }
                            else {
                                roff[v][u][t][s] = roff[v - 1][u][t][s];
                                free(ptr);
                            }
                        }
                    }
                }
            }

            break;
        default:
            fprintf(stderr, "Wrong usage of init_ref_offset function!\n");
            exit(-10);
    }

    return (roff);
}

void free_ref_offset(int vu[2], int ts[2], int type, int prd_order, int *****roff) {
    int dt, ds, v, u, t, s, k;

    int min_ds, max_ds, min_dt, max_dt;
    int mmin_ds, mmax_ds, mmin_dt, mmax_dt;

    min_ds = max_ds = min_dt = max_dt = 0;
    mmin_ds = mmax_ds = mmin_dt = mmax_dt = 0;

    switch (type) {
        case INTRA_PRED: // Intra
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
            for (v = 0; v < vu[HEIGHT]; v++) {
                for (u = 0; u < vu[WIDTH]; u++) {
                    for (t = 0; t < ts[HEIGHT]; t++) {
                        for (s = 0; s < ts[WIDTH]; s++) {
                            //Conditions to check which references are available for each pixel
                            if (((v == 0 && u == 0) && ((t == 0 && (s == 0 || s + min_ds <= 0 || s + max_ds >= ts[WIDTH])) ||
                                                        (t + min_dt <= 0 && (s == 0 || s + min_ds <= 0 || s + max_ds >= ts[WIDTH])))) ||
                                ((v > 0 || u > 0) && (t == 0 && s == 0))) {
                                free(roff[v][u][t][s]);
                            }
                        }
                    }
                }
            }

            safefree((void **) &(roff));

            break;

        case MI_UP_PRED:
            for (k = 0; k < prd_order; k++) {
                dt = dyx[k].y;
                ds = dyx[k].x;

                if (dt < min_dt) min_dt = dt;
                if (ds < min_ds) min_ds = ds;
                if (ds > max_ds) max_ds = ds;
            }

            for (k = 0; k < prd_order - 1; k++) {
                dt = idyx[k].y;
                ds = idyx[k].x;

                if (dt < mmin_dt) mmin_dt = dt;
                if (dt > mmax_dt) mmax_dt = dt;
                if (ds < mmin_ds) mmin_ds = ds;
                if (ds > mmax_ds) mmax_ds = ds;
            }

            // Cycle that runs for all the pixels
            for (v = 0; v < vu[HEIGHT]; v++) {
                for (u = 0; u < vu[WIDTH]; u++) {
                    for (t = 0; t < ts[HEIGHT]; t++) {
                        for (s = 0; s < ts[WIDTH]; s++) {
                            if (v == 0 && u == 0) {
                                if ((t == 0 && (s == 0 || s + min_ds <= 0 || s + max_ds >= ts[WIDTH])) ||
                                    (t + min_dt <= 0 && (s == 0 || s + min_ds <= 0 || s + max_ds >= ts[WIDTH]))) {
                                    free(roff[v][u][t][s]);
                                }
                            }
                            else if (v == 0 && t == 0 && s == 0) {
                                free(roff[v][u][t][s]);
                            }
                            else if (v == 1 && u == 0) {
                                //Conditions to check which references are available for each pixel
                                if ((t == 0 && (s == 0 || s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH])) ||
                                    (t + mmin_dt <= 0 && (s == 0 || s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH])) ||
                                    (t + mmax_dt >= ts[HEIGHT] && (s == 0 || s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH]))) {
                                    free(roff[v][u][t][s]);
                                }
                            }
                        }
                    }
                }
            }

            safefree((void **) &(roff));

            break;

        case MI_LEFT_PRED:
            for (k = 0; k < prd_order; k++) {
                dt = dyx[k].y;
                ds = dyx[k].x;

                if (dt < min_dt) min_dt = dt;
                if (ds < min_ds) min_ds = ds;
                if (ds > max_ds) max_ds = ds;
            }

            for (k = 0; k < prd_order - 1; k++) {
                dt = tridyx[k].y;
                ds = tridyx[k].x;

                if (dt < mmin_dt) mmin_dt = dt;
                if (dt > mmax_dt) mmax_dt = dt;
                if (ds < mmin_ds) mmin_ds = ds;
                if (ds > mmax_ds) mmax_ds = ds;
            }

            // Cycle that runs for all the pixels
            for (v = 0; v < vu[HEIGHT]; v++) {
                for (u = 0; u < vu[WIDTH]; u++) {
                    for (t = 0; t < ts[HEIGHT]; t++) {
                        for (s = 0; s < ts[WIDTH]; s++) {
                            if (v == 0 && u == 0) {
                                if ((t == 0 && (s == 0 || s + min_ds <= 0 || s + max_ds >= ts[WIDTH])) ||
                                    (t + min_dt <= 0 && (s == 0 || s + min_ds <= 0 || s + max_ds >= ts[WIDTH]))) {
                                    free(roff[v][u][t][s]);
                                }
                            }
                            else if (u == 0 && t == 0 && s == 0) {
                                free(roff[v][u][t][s]);
                            }
                            else if (v == 0 && u == 1) {
                                //Conditions to check which references are available for each pixel
                                if ((t == 0 && (s == 0 || s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH])) ||
                                    (t + mmin_dt <= 0 && (s == 0 || s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH]))) {
                                    free(roff[v][u][t][s]);
                                }
                            }
                        }
                    }
                }
            }

            safefree((void **) &(roff));

            break;

        case MI_LDIAG_PRED:
            for (k = 0; k < prd_order; k++) {
                dt = dyx[k].y;
                ds = dyx[k].x;

                if (dt < min_dt) min_dt = dt;
                if (ds < min_ds) min_ds = ds;
                if (ds > max_ds) max_ds = ds;
            }

            for (k = 0; k < prd_order - 1; k++) {
                dt = idyx[k].y;
                ds = idyx[k].x;

                if (dt < mmin_dt) mmin_dt = dt;
                if (dt > mmax_dt) mmax_dt = dt;
                if (ds < mmin_ds) mmin_ds = ds;
                if (ds > mmax_ds) mmax_ds = ds;
            }

            // Cycle that runs for all the pixels
            for (v = 0; v < vu[HEIGHT]; v++) {
                for (u = 0; u < vu[WIDTH]; u++) {
                    for (t = 0; t < ts[HEIGHT]; t++) {
                        for (s = 0; s < ts[WIDTH]; s++) {
                            if (v == 0 || u == 0) {
                                if ((t == 0 && (s == 0 || s + min_ds <= 0 || s + max_ds >= ts[WIDTH])) ||
                                    (t + min_dt <= 0 && (s == 0 || s + min_ds <= 0 || s + max_ds >= ts[WIDTH]))) {
                                    free(roff[v][u][t][s]);
                                }
                            }
                            else if ((v == 0 || u == 0) && t == 0 && s == 0) {
                                free(roff[v][u][t][s]);
                            }
                            else if (v == 0) {
                                if (t == 0 && s == 0) {

                                }
                            }
                            else if (v == 1 && u == 1) {
                                //Conditions to check which references are available for each pixel
                                if ((t == 0 && (s == 0 || s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH])) ||
                                    (t + mmin_dt <= 0 && ( s == 0 || s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH])) ||
                                    (t + mmax_dt >= ts[HEIGHT] && (s == 0 || s + mmin_ds <= 0 || s + max_ds >= ts[WIDTH]))) {
                                    free(roff[v][u][t][s]);
                                }
                            }
                        }
                    }
                }
            }
            safefree((void **) &(roff));

            break;

        case MI_RDIAG_PRED:
            for (k = 0; k < prd_order; k++) {
                dt = dyx[k].y;
                ds = dyx[k].x;

                if (dt < min_dt) min_dt = dt;
                if (ds < min_ds) min_ds = ds;
                if (ds > max_ds) max_ds = ds;
            }

            for (k = 0; k < prd_order - 1; k++) {
                dt = idyx[k].y;
                ds = idyx[k].x;

                if (dt < mmin_dt) mmin_dt = dt;
                if (dt > mmax_dt) mmax_dt = dt;
                if (ds < mmin_ds) mmin_ds = ds;
                if (ds > mmax_ds) mmax_ds = ds;
            }

            // Cycle that runs for all the pixels
            for (v = 0; v < vu[HEIGHT]; v++) {
                for (u = 0; u < vu[WIDTH]; u++) {
                    for (t = 0; t < ts[HEIGHT]; t++) {
                        for (s = 0; s < ts[WIDTH]; s++) {
                            if (v == 0 || u == vu[WIDTH] - 1) {
                                if ((t == 0 && (s == 0 || s + min_ds <= 0 || s + max_ds >= ts[WIDTH])) ||
                                    (t + min_dt <= 0 && (s == 0 || s + min_ds <= 0 || s + max_ds >= ts[WIDTH]))) {
                                    free(roff[v][u][t][s]);
                                }
                            }
                            else if (v == 0 && t == 0 && s == 0) {
                                free(roff[v][u][t][s]);
                            }
                            else if (v == 1 && u == 0) {
                                //Conditions to check which references are available for each pixel
                                if ((t == 0 && (s == 0 || s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH])) ||
                                    (t + mmin_dt <= 0 && ( s == 0 || s + mmin_ds <= 0 || s + mmax_ds >= ts[WIDTH])) ||
                                    (t + mmax_dt >= ts[HEIGHT] && (s == 0 || s + mmin_ds <= 0 || s + max_ds >= ts[WIDTH]))) {
                                    free(roff[v][u][t][s]);
                                }
                            }
                        }
                    }
                }
            }
            safefree((void **) &(roff));

            break;

        default:
            fprintf(stderr, "Wrong usage of free_ref_offset function!\n");
            exit(-11);
    }
}

int *init_ctx_weight(int type, int prd_order, int delta) {
    int *ctx_weight, k;
    double dy, dx;

    ctx_weight = (int *) alloc_mem((prd_order) * sizeof(int));

    switch (type) {
        case INTRA_PRED: // Intra
        case MI_BORDERS_PRED:
            for (k = 0; k < prd_order; k++) {
                dy = dyx[k].y;
                dx = dyx[k].x;

                ctx_weight[k] = (int) (64.0 / sqrt(dy * dy + dx * dx) + 0.5);
            }

            break;

        case MI_UP_PRED:
        case MI_LDIAG_PRED:
        case MI_RDIAG_PRED:
            ctx_weight[0] = (int) (64.0 / sqrt(delta * delta) + 0.5);

            for (k = 0; k < prd_order - 1; k++) {
                dy = idyx[k].y;
                dx = idyx[k].x;

                ctx_weight[k + 1] = (int) (64.0 / sqrt(delta * delta + dy * dy + dx * dx) + 0.5);
            }

            break;

        case MI_LEFT_PRED:
            ctx_weight[0] = (int) (64.0 / sqrt(delta * delta) + 0.5);

            for (k = 0; k < prd_order - 1; k++) {
                dy = tridyx[k].y;
                dx = tridyx[k].x;

                ctx_weight[k + 1] = (int) (64.0 / sqrt(delta * delta + dy * dy + dx * dx) + 0.5);
            }

            break;

        default:
            fprintf(stderr, "Wrong usage of init_ctx_weight function!\n");
            exit(-12);
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

void mtf_classlabel(char ****class, int *mtfbuf, int v, int u, int t, int s, int bsize, int width, int num_class) {
    int i, j, k, ref[5];

    if (v == 0) {
        if (u == 0) {
            if (t == 0) {
                if (s == 0) {
                    ref[0] = ref[1] = ref[2] = ref[3] = ref[4] = 0;
                }
                else {
                    ref[0] = ref[1] = ref[2] = ref[3] = ref[4] = class[v][u][t][s - 1];
                }
            }
            else {
                ref[0] = ref[3] = class[v][u][t - 1][s];
                ref[1] = ref[4] = (s == 0) ? class[v][u][t - 1][s] : class[v][u][t][s - 1];
                ref[2] = (s + bsize >= width) ? class[v][u][t - 1][s] : class[v][u][t - 1][s + bsize];

                if (ref[1] == ref[2]) {
                    ref[2] = ref[0];
                    ref[0] = ref[1];
                }
            }
        }
        else {
            ref[0] = (t == 0 && s == 0) ? class[v][u - 1][t][s] : (t == 0) ? class[v][u][t][s - 1] : class[v][u][t - 1][s];
            ref[1] = (t == 0 && s == 0) ? class[v][u - 1][t][s] : (s == 0) ? class[v][u][t - 1][s] : class[v][u][t][s - 1];
            ref[2] = (t == 0 || s + bsize >= width) ? class[v][u - 1][t][s] : class[v][u][t - 1][s + bsize];
            ref[3] = ref[4] = class[v][u - 1][t][s];

            if (ref[1] == ref[2]) {
                ref[2] = ref[0];
                ref[0] = ref[1];
            }
        }
    }
    else {
        ref[0] = (t == 0 && s == 0) ? class[v - 1][u][t][s] : (t == 0) ? class[v][u][t][s - 1] : class[v][u][t - 1][s];
        ref[1] = (t == 0 && s == 0) ? class[v - 1][u][t][s] : (s == 0) ? class[v][u][t - 1][s] : class[v][u][t][s - 1];
        ref[2] = (t == 0 || s + bsize >= width) ? class[v - 1][u][t][s] : class[v][u][t - 1][s + bsize];
        ref[3] = (u == 0) ? class[v - 1][u][t][s] : class[v][u - 1][t][s];
        ref[4] = class[v - 1][u][t][s];

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

int **copy_bref(int bframes, int bref_const[][5]) {
    int i, j;

    int **bref = (int **) alloc_2d_array(bframes, 5, sizeof(int));

    for (i = 0; i < bframes; i++) {
        for (j = 0; j < 5; j++) {
            bref[i][j] = bref_const[i][j];
        }
    }

    return bref;
}

int **select_bref(int bframes) {
    int **bref = NULL;

    switch(bframes) {
        case 2:
            bref = copy_bref(bframes, bref1);
            break;
        case 3:
            bref = copy_bref(bframes, bref2);
            break;
        case 4:
            bref = copy_bref(bframes, bref3);
            break;
        case 5:
            bref = copy_bref(bframes, bref4);
            break;
        case 6:
            bref = copy_bref(bframes, bref5);
            break;
        case 7:
            bref = copy_bref(bframes, bref6);
            break;
        case 8:
            bref = copy_bref(bframes, bref7);
            break;
        case 9:
            bref = copy_bref(bframes, bref8);
            break;
        case 10:
            bref = copy_bref(bframes, bref9);
            break;
        default:
            bref = copy_bref(bframes, bref3);
            break;
    }

    return bref;
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
