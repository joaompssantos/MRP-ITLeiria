/*
“Commons Clause” License Condition v1.0

The Software is provided to you by the Licensor under the License, as defined below, subject to the following condition.

Without limiting other conditions in the License, the grant of rights under the License will not include, and the License does not grant to you, the right to Sell the Software.

For purposes of the foregoing, “Sell” means practicing any or all of the rights granted to you under the License to provide to third parties, for a fee or other consideration (including without limitation fees for hosting or consulting/ support services related to the Software), a product or service whose value derives, entirely or substantially, from the functionality of the Software. Any license notice or attribution required by the License must also include this Commons Clause License Condition notice.

Software: Minimum Rate Predictors

License: BSD-3-Clause with Commons Clause

Licensor: Instituto de Telecomunicações

Copyright © 2023 Instituto de Telecomunicações. All Rights Reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

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

const int mask[9] = {0x00FF,
				 	 0x01FF,
					 0x03FF,
					 0x07FF,
					 0x0FFF,
					 0x1FFF,
					 0x3FFF,
					 0x7FFF,
					 0xFFFF
};

double sigma_h[] = {0.85, 1.15, 1.50, 1.90, 2.55, 3.30, 4.25, 5.60, 7.15, 9.20,12.05,15.35,19.95,25.85,32.95,44.05};
double sigma_a[] = {0.15, 0.26, 0.38, 0.57, 0.83, 1.18, 1.65, 2.31, 3.22, 4.47, 6.19, 8.55,11.80,16.27,22.42,30.89};
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

void safefree_yuv(IMAGE **pp) {
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
        fprintf(stderr, "Can\'t allocate memory (size = %d)!\n", (int) size);
        exit(1);
    }
    return (ptr);
}

// Image matrix allocation
void **alloc_2d_array(int height, int width, int size) {
    void **mat;
    char *ptr;
    int k;

    mat = (void **) alloc_mem(sizeof(void **) * height + height * width * size);
    ptr = (char *) (mat + height);

    for (k = 0; k < height; k++) {
        mat[k] = ptr;
        ptr += width * size;
    }

    return (mat);
}

// Image matrix allocation
void ***alloc_3d_array(int height, int width, int frames, int size) {
    void ***mat;
    int i, j;

    mat = (void ***) alloc_mem(
            sizeof(void ***) * frames + sizeof(void **) * frames * height + size * frames * height * width);

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

// Function to alloc image
IMAGE *alloc_image(int width, int height, int maxval) {
    // Image struct vector
    IMAGE *img;

    // Image struct memory allocation
    img = (IMAGE *) alloc_mem(sizeof(IMAGE));

    // Set image values
    img->width = width;
    img->height = height;
    img->maxval = maxval;

    img->val = (img_t **) alloc_2d_array(img->height, img->width, sizeof(img_t));
    return (img);
}

// Function to copy YUV image
IMAGE *copy_yuv(IMAGE *img) {
    int i, j;

    IMAGE *new_img = alloc_image(img->width, img->height, img->maxval);

    for (i = 0; i < new_img->height; i++) {
        for (j = 0; j < new_img->width; j++) {
            new_img->val[i][j] = img->val[i][j];
        }
    }

    return (new_img);
}

// Function to reverse the endianness of a unsigned short
unsigned short reverse_endianness(unsigned short s, int endianness) {
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
void write_yuv(IMAGE *img, char *filename, int depth, int endianness) {
    int i, j;
    unsigned short byte;
    FILE *fp;

    fp = fileopen(filename, "ab");

    for (i = 0; i < img->height; i++) {
        for (j = 0; j < img->width; j++) {
            if (depth == 8) {
                putc(img->val[i][j], fp);
            }
            else if (depth > 8) {
                byte = reverse_endianness(img->val[i][j], endianness);

                putc((byte >> 8) & 0x00FF, fp);
                putc(byte & 0x00FF, fp);
            }
        }
    }

/*	for (i = 0; i < img->height / 2; i++) {
		for (j = 0; j < img->width / 2; j++) {
			if (depth > 8) {
				byte = reverse_endianness((int) (pow(2, depth) / 2), endianness);

				putc((byte >> 8) & 0x00FF, fp);
				putc(byte & 0x00FF, fp);
			}
			else {
				putc((int) (pow(2, depth) / 2), fp);
			}
		}
	}

	for (i = 0; i < img->height / 2; i++) {
		for (j = 0; j < img->width / 2; j++) {
			if (depth > 8) {
				byte = reverse_endianness((int) (pow(2, depth) / 2), endianness);

				putc((byte >> 8) & 0x00FF, fp);
				putc(byte & 0x00FF, fp);
			}
			else {
				putc((int) (pow(2, depth) / 2), fp);
			}
		}
	}*/

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
 |
 |  Returns:  IMAGE* --> returns a video type structure
 *-------------------------------------------------------------------*/
IMAGE *read_yuv(char *filename, int height, int width, int frame, int depth, int endianness, int chroma) {
    int i, j, shift_first, shift_second, first, second;
    double chroma_pass = (chroma == GRAY ? 1 : (chroma == S444 ? 3 : 1.5));
    IMAGE *img;
    FILE *fp;

    //Open file
    fp = fileopen(filename, "rb");

    // Check if image dimensions are correct (It has to be multiple of BASE_BSIZE)
    if ((width % BASE_BSIZE) || (height % BASE_BSIZE)) {
        fprintf(stderr, "Image width and height must be multiples of %d!\n", BASE_BSIZE);
        exit(1);
    }

    if (endianness == LITTLE_ENDIANNESS) {
        shift_first = 0;
        shift_second = 8;
    }
    else {
        shift_first = 8;
        shift_second = 0;
    }

    // Image allocation
    img = alloc_image(width, height, (int) (pow(2, depth) - 1));

    if (frame > 0) {
        if (depth > 8) {
            fseek(fp, (long) (img->height * img->width * 2 * chroma_pass * frame), SEEK_SET);
        }
        else {
            fseek(fp, (long) (img->height * img->width * chroma_pass * frame), SEEK_SET);
        }
    }

    for (i = 0; i < img->height; i++) {
        for (j = 0; j < img->width; j++) {
            first = (img_t) fgetc(fp);

            if (depth > 8) {
                second = (img_t) fgetc(fp);

                img->val[i][j] = (img_t) ((first << shift_first) + (second << shift_second));
                img->val[i][j] = (img_t) (img->val[i][j] & mask[depth - 8]);
            }
            else {
                img->val[i][j] = (img_t) first;
            }
        }
    }

/*	for (i = 0; i < img->height / 2; i++) {
		for (j = 0; j < img->width / 2; j++) {
			fgetc(fp);

			if (depth > 8) {
				fgetc(fp);
			}
		}
	}

	for (i = 0; i < img->height / 2; i++) {
		for (j = 0; j < img->width / 2; j++) {
			fgetc(fp);

			if (depth > 8) {
				fgetc(fp);
			}
		}
	}*/

    fclose(fp);
    return (img);
}

/*--------------------------- init_ref_offset -----------------------*
 |  Function init_ref_offset
 |
 |  Purpose:  Calculates the offset of a reference to a given pixel
 |            and a specific type
 |
 |  Parameters:
 |		height          --> Number of height pixels (IN)
 |		width           --> Number of width pixels (IN)
 |		type            --> Type of offset to calculate (IN)
 |		prd_order	    --> Number of reference pixels for type (IN)
 |		mi_size     	--> Micro image size (IN)
 |
 |  Returns:  int*** --> Array with the references offset
 *-------------------------------------------------------------------*/
int ***init_ref_offset(int height, int width, int type, int prd_order) {
    int ***roff, *ptr;
    int x, y, dx, dy, k, base = 0;

    int min_dx, max_dx, min_dy, max_dy;

    min_dx = max_dx = min_dy = max_dy = 0;

    roff = (int ***) alloc_2d_array(height, width, sizeof(int *));

    switch (type) {
        case INTRA_PRED: // Intra
            base = 0;
            break;

        case BACK_PRED: // Inter backwards
            base = -width * (height + 1);
            break;

        case FOR_PRED: // Inter forward
            base = width * (height + 1);
            break;

        default:
            fprintf(stderr, "Wrong usage of init_ref_offset function!\n");
            exit(-10);
    }

    switch (type) {
        case INTRA_PRED: // Intra
            //Values to check for special cases
            for (k = 0; k < prd_order; k++) {
                dy = dyx[k].y;
                dx = dyx[k].x;

                if (dy < min_dy) min_dy = dy;
                if (dx < min_dx) min_dx = dx;
                if (dy > max_dy) max_dy = dy;
                if (dx > max_dx) max_dx = dx;
            }

            //Cycle that runs for all the pixels
            for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                    ptr = (int *) alloc_mem((prd_order) * sizeof(int));
                    //Conditions to check which references are available for each pixel
                    if (y == 0) {
                        if (x == 0) {
                            roff[y][x] = ptr;
                            dx = 0;
                            dy = height;

                            for (k = 0; k < prd_order; k++) {
                                *ptr++ = dy * width + dx; //Points to a line filled with 128 (if max_val = 256)
                            }
                        }
                        else if (x + min_dx <= 0 || x + max_dx >= width) {
                            roff[y][x] = ptr;
                            dy = 0;

                            for (k = 0; k < prd_order; k++) {
                                dx = dyx[k].x;

                                if (x + dx < 0) dx = -x;
                                else if (dx >= 0) dx = -1;

                                *ptr++ = dy * width + dx;
                            }
                        }
                        else {
                            roff[y][x] = roff[y][x - 1];
                            free(ptr);
                        }
                    }
                    else if (y + min_dy <= 0) {
                        if (x == 0) {
                            roff[y][x] = ptr;

                            for (k = 0; k < prd_order; k++) {
                                dy = dyx[k].y;

                                if (y + dy < 0) dy = -y;
                                else if (dy >= 0) dy = -1;

                                dx = dyx[k].x;

                                if (x + dx < 0) dx = -x;

                                *ptr++ = dy * width + dx;
                            }
                        }
                        else if (x + min_dx <= 0 || x + max_dx >= width) {
                            roff[y][x] = ptr;

                            for (k = 0; k < prd_order; k++) {
                                dy = dyx[k].y;

                                if (y + dy < 0) dy = -y;
                                dx = dyx[k].x;

                                if (x + dx < 0) dx = -x;
                                else if (x + dx >= width) {
                                    dx = width - x - 1;
                                }

                                *ptr++ = dy * width + dx;
                            }
                        }
                        else {
                            roff[y][x] = roff[y][x - 1];
                            free(ptr);
                        }
                    }
                    else {
                        roff[y][x] = roff[y - 1][x];
                        free(ptr);
                    }
                }
            }

            break;

        case BACK_PRED: // Inter Backwards
        case FOR_PRED: // Inter Forward
            for (k = 0; k < prd_order - 1; k++) {
                dy = idyx[k].y;
                dx = idyx[k].x;

                if (dy < min_dy) min_dy = dy;
                if (dy > max_dy) max_dy = dy;
                if (dx < min_dx) min_dx = dx;
                if (dx > max_dx) max_dx = dx;
            }

            //Cycle that runs for all the pixels
            for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                    ptr = (int *) alloc_mem((prd_order) * sizeof(int));

                    if (y == 0) {
                        if (x == 0) {
                            roff[y][x] = ptr;
                            *ptr++ = base;

                            for (k = 0; k < prd_order - 1; k++) {
                                dy = idyx[k].y;
                                dx = idyx[k].x;

                                if (y + dy < 0 || x + dx < 0) {
                                    *ptr++ = base;
                                }
                                else {
                                    *ptr++ = dy * width + dx + base;
                                }
                            }
                        }
                        else if (x + min_dx <= 0 || x + max_dx >= width) {
                            roff[y][x] = ptr;
                            *ptr++ = base;

                            for (k = 0; k < prd_order - 1; k++) {
                                dy = idyx[k].y;
                                dx = idyx[k].x;

                                if (y + dy < 0 || x + dx < 0 || x + dx >= width) {
                                    *ptr++ = base;
                                }
                                else {
                                    *ptr++ = dy * width + dx + base;
                                }
                            }
                        }
                        else {
                            roff[y][x] = roff[y][x - 1];
                            free(ptr);
                        }
                    }
                    else if (y + min_dy <= 0) {
                        if (x == 0) {
                            roff[y][x] = ptr;
                            *ptr++ = base;

                            for (k = 0; k < prd_order - 1; k++) {
                                dy = idyx[k].y;
                                dx = idyx[k].x;

                                if (y + dy < 0 || x + dx < 0) {
                                    *ptr++ = base;
                                }
                                else {
                                    *ptr++ = dy * width + dx + base;
                                }
                            }
                        }
                        else if (x + min_dx <= 0 || x + max_dx >= width) {
                            roff[y][x] = ptr;
                            *ptr++ = base;

                            for (k = 0; k < prd_order - 1; k++) {
                                dy = idyx[k].y;
                                dx = idyx[k].x;

                                if (y + dy < 0 || x + dx < 0 || x + dx >= width) {
                                    *ptr++ = base;
                                }
                                else {
                                    *ptr++ = dy * width + dx + base;
                                }
                            }
                        }
                        else {
                            roff[y][x] = roff[y][x - 1];
                            free(ptr);
                        }
                    }
                    else if (y + max_dy >= height) {
                        if (x == 0 || x + min_dx <= 0 || x + max_dx >= width) {
                            roff[y][x] = ptr;
                            *ptr++ = base;

                            for (k = 0; k < prd_order - 1; k++) {
                                dy = idyx[k].y;
                                dx = idyx[k].x;

                                if (y + dy >= height || x + dx < 0 || x + dx >= width) {
                                    *ptr++ = base;
                                }
                                else {
                                    *ptr++ = dy * width + dx + base;
                                }
                            }
                        }
                        else {
                            roff[y][x] = roff[y][x - 1];
                            free(ptr);
                        }
                    }
                    else {
                        roff[y][x] = roff[y - 1][x];
                        free(ptr);
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

/*--------------------------- free_ref_offset -----------------------*
 |  Function free_ref_offset
 |
 |  Purpose:  Free the offset reference array
 |
 |  Parameters:
 |		height          --> Number of height pixels (IN)
 |		width           --> Number of width pixels (IN)
 |		type            --> Type of offset (IN)
 |		prd_order	    --> Number of reference pixels for type (IN)
 |		roff            --> Refernce offset array to free (IN)
 |
 |  Returns:  void
 *-------------------------------------------------------------------*/
void free_ref_offset(int height, int width, int type, int prd_order, int ***roff) {
    int dy, dx, x, y, k;

    int min_dx, max_dx, min_dy, max_dy;

    min_dx = max_dx = min_dy = max_dy = 0;

    switch (type) {
        case INTRA_PRED: // Intra
            // Check if mi_size should be also used
            for (k = 0; k < prd_order; k++) {
                dy = dyx[k].y;
                dx = dyx[k].x;

                if (dy < min_dy) min_dy = dy;
                if (dy > max_dy) max_dy = dy;
                if (dx < min_dx) min_dx = dx;
                if (dx > max_dx) max_dx = dx;
            }

            //Cycle that runs for all the pixels
            for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                    //Conditions to check which references are available for each pixel
                    if ((y == 0 && (x == 0 || x + min_dx <= 0 || x + max_dx >= width)) ||
                        (y + min_dy <= 0 && (x == 0 || x + min_dx <= 0 || x + max_dx >= width))) {
                        free(roff[y][x]);
                    }
                }
            }

            safefree((void **) &(roff));

            break;

        case BACK_PRED: // Inter Backwards
        case FOR_PRED: // Inter Forward
            // Check if mi_size should be also used
            for (k = 0; k < prd_order - 1; k++) {
                dy = idyx[k].y;
                dx = idyx[k].x;

                if (dy < min_dy) min_dy = dy;
                if (dy > max_dy) max_dy = dy;
                if (dx < min_dx) min_dx = dx;
                if (dx > max_dx) max_dx = dx;
            }

            //Cycle that runs for all the pixels
            for (y = 0; y < height; y++) {
                for (x = 0; x < width; x++) {
                    //Conditions to check which references are available for each pixel
                    if ((y == 0 && (x == 0 || x + min_dx <= 0 || x + max_dx >= width)) ||
                        (y + min_dy <= 0 && (x == 0 || x + min_dx <= 0 || x + max_dx >= width)) ||
                        (y + max_dy >= height && (x == 0 || x + min_dx <= 0 || x + max_dx >= width))) {
                        free(roff[y][x]);
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

/*--------------------------- init_ctx_weight -----------------------*
 |  Function init_ctx_weight
 |
 |  Purpose:  Calculates the ctx weights for each reference pixel of
 |            of a given type
 |
 |  Parameters:
 |		type            --> Type of pred to calculate (IN)
 |		prd_order	    --> Number of reference pixels for type (IN)
 |		delta     	    --> Micro image size (IN)
 |
 |  Returns:  int*      --> Array with the ctx weights
 *-------------------------------------------------------------------*/
int *init_ctx_weight(int type, int prd_order, int delta) {
    int *ctx_weight, k;
    double dy, dx;

    ctx_weight = (int *) alloc_mem((prd_order) * sizeof(int));

    switch (type) {
        case INTRA_PRED: // Intra
            for (k = 0; k < prd_order; k++) {
                dy = dyx[k].y;
                dx = dyx[k].x;

                ctx_weight[k] = (int) (64.0 / sqrt(dy * dy + dx * dx) + 0.5);
            }

            break;

        case BACK_PRED: // Inter backwards
        case FOR_PRED: // Inter forward
            ctx_weight[0] = (int) (64.0 / sqrt(delta * delta) + 0.5);

            for (k = 0; k < prd_order - 1; k++) {
                dy = idyx[k].y;
                dx = idyx[k].x;

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

void mtf_classlabel(char **class, int *mtfbuf, int y, int x, int bsize, int width, int num_class) {
    int i, j, k, ref[3];

    if (y == 0) {
        if (x == 0) {
            ref[0] = ref[1] = ref[2] = 0;
        }
        else {
            ref[0] = ref[1] = ref[2] = class[y][x - 1];
        }
    }
    else {
        ref[0] = class[y - 1][x];
        ref[1] = (x == 0) ? class[y - 1][x] : class[y][x - 1];
        ref[2] = (x + bsize >= width) ? class[y - 1][x] : class[y - 1][x + bsize];
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

    switch (bframes) {
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

char *cat_str(char *str1, char *str2, int type) {
    char *new_str = (char *) alloc_mem(sizeof(char) * (strlen(str1) + strlen(str2) + 1));

    if (new_str != NULL) {
        new_str[0] = '\0';   // ensures the memory is an empty string

        strcat(new_str, str1);
        strcat(new_str, str2);
    }
    else {
        printf("malloc failed!\n");
        exit(-12);
    }

    if (type == 1) safefree((void **) &str1);
    safefree((void **) &str2);

    return new_str;
}

// buffer must have length >= sizeof(int) + 1
// Write to the buffer backwards so that the binary representation
// is in the correct order i.e.  the LSB is on the far right
// instead of the far left of the printed string
char *int2bin(int a, int buf_size) {
    char *buffer = (char *) alloc_mem(sizeof(char) * (buf_size + 1));
    buffer[buf_size] = '\0';

    for (int i = buf_size - 1; i >= 0; i--) {
        buffer[i] = (char) ((a & 1) + '0');

        a >>= 1;
    }

    return buffer;
}
