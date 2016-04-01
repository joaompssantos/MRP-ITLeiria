#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mrp.h"

//Reference pixel position separated by their distances
const POINT dyx[] = {
		/* 1 */
		{ 0,-1}, {-1, 0},
		/* 2 */
		{ 0,-2}, {-1,-1}, {-2, 0}, {-1, 1},
		/* 3 */
		{ 0,-3}, {-1,-2}, {-2,-1}, {-3, 0}, {-2, 1}, {-1, 2},
		/* 4 */
		{ 0,-4}, {-1,-3}, {-2,-2}, {-3,-1}, {-4, 0}, {-3, 1}, {-2, 2}, {-1, 3},
		/* 5 */
		{ 0,-5}, {-1,-4}, {-2,-3}, {-3,-2}, {-4,-1}, {-5, 0}, {-4, 1}, {-3, 2}, {-2, 3}, {-1, 4},
		/* 6 */
		{ 0,-6}, {-1,-5}, {-2,-4}, {-3,-3}, {-4,-2}, {-5,-1}, {-6, 0}, {-5, 1}, {-4, 2}, {-3, 3}, {-2, 4}, {-1, 5},
		/* 7 */
		{ 0,-7}, {-1,-6}, {-2,-5}, {-3,-4}, {-4,-3}, {-5,-2}, {-6,-1}, {-7, 0}, {-6, 1}, {-5, 2}, {-4, 3}, {-3, 4}, {-2, 5}, {-1, 6},
		/* 8 */
		{ 0,-8}, {-1,-7}, {-2,-6}, {-3,-5}, {-4,-4}, {-5,-3}, {-6,-2}, {-7,-1}, {-8, 0}, {-7, 1}, {-6, 2}, {-5, 3}, {-4, 4}, {-3, 5}, {-2, 6}, {-1, 7},
};
//Reference pixel position separated by their distances
// {y, x}
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
//		/* 5 */
//		{ 0,-5}, {-1,-4}, {-2,-3}, {-3,-2}, {-4,-1}, {-5, 0}, {-4, 1}, {-3, 2},
//		{-2, 3}, {-1, 4},
//		/* 6 */
//		{ 0,-6}, {-1,-5}, {-2,-4}, {-3,-3}, {-4,-2}, {-5,-1}, {-6, 0}, {-5, 1},
//		{-4, 2}, {-3, 3}, {-2, 4}, {-1, 5},
//		/* 7 */
//		{ 0,-7}, {-1,-6}, {-2,-5}, {-3,-4}, {-4,-3}, {-5,-2}, {-6,-1}, {-7, 0},
//		{-6, 1}, {-5, 2}, {-4, 3}, {-3, 4}, {-2, 5}, {-1, 6},
//		/* 8 */
//		{ 0,-8}, {-1,-7}, {-2,-6}, {-3,-5}, {-4,-4}, {-5,-3}, {-6,-2}, {-7,-1},
//		{-8, 0}, {-7, 1}, {-6, 2}, {-5, 3}, {-4, 4}, {-3, 5}, {-2, 6}, {-1, 7},
};

int bref1[4][5] = {{2, -2, -1, 0, -1},
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

const mask[9] = {0x00FF,
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
		fprintf(stderr, "Can\'t allocate memory (size = %d)!\n", (int)size);
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
		mat[k] =  ptr;
		ptr += width * size;
	}

	return (mat);
}

// Image matrix allocation
void ***alloc_3d_array(int height, int width, int frames, int size) {
	void ***mat;
	int i, j;

	mat = (void ***) alloc_mem(sizeof(void ***) * frames + sizeof(void **) * frames * height + size * frames  * height * width);

	void** hei = (void**)(mat + frames);
	char* wid = (char*)(hei + frames * height);

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
	img = (IMAGE *)alloc_mem(sizeof(IMAGE));

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

	return(new_img);
}

// Function to reverse the endianness of a unsigned short
unsigned short reverse_endianness (unsigned short s, int endianness) {
    unsigned char c1, c2;

    if (endianness == 1) {
        return s;
    }
    else {
        c1 = s & 255;
        c2 = (s >> 8) & 255;

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

	for (i = 0; i < img->height / 2; i++) {
		for (j = 0; j < img->width / 2; j++) {
			putc((int) (pow(2, depth) / 2), fp);

			if (depth > 8) {
				putc((int) (pow(2, depth) / 2), fp);
			}
		}
	}

	for (i = 0; i < img->height / 2; i++) {
		for (j = 0; j < img->width / 2; j++) {
			putc((int) (pow(2, depth) / 2), fp);

			if (depth > 8) {
				putc((int) (pow(2, depth) / 2), fp);
			}
		}
	}

	fclose(fp);

	return;
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
 |
 |  Returns:  IMAGE* --> returns a video type structure
 *-------------------------------------------------------------------*/
IMAGE *read_yuv(char *filename, int height, int width, int frame, int depth, int endianness) {
	int i, j, shift_first, shift_second, first, second;
	IMAGE *img;
	FILE *fp;

	//Open file
	fp = fileopen(filename, "rb");

	// Check if image dimensions are correct (It has to be multiple of BASE_BSIZE)
	if ((width % BASE_BSIZE) || (height % BASE_BSIZE)) {
		fprintf(stderr, "Image width and height must be multiples of %d!\n", BASE_BSIZE);
		exit(1);
	}

	if( endianness == LITTLE_ENDIANNESS ) {
		shift_first = 0;
		shift_second = 8;
	}
	else if (endianness == BIG_ENDIANNESS){
		shift_first = 8;
		shift_second = 0;
	}

	// Image allocation
	img = alloc_image(width, height, (int) (pow(2, depth) - 1));

	if (frame > 0) fseek(fp, img->height * img->width * 1.5 * frame, SEEK_SET);

	for (i = 0; i < img->height; i++) {
		for (j = 0; j < img->width; j++) {
			first = (img_t)fgetc(fp);

			if (depth > 8) {
				second = (img_t)fgetc(fp);

				img->val[i][j] = (first << shift_first) + (second << shift_second);
				img->val[i][j] = img->val[i][j] & mask[depth - 8];
			}
//			printf("Testes: %d\n", img->val[i][j]);getchar();
		}
	}

	for (i = 0; i < img->height / 2; i++) {
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
	}

	fclose(fp);
	return (img);
}


int *gen_hufflen(uint *hist, int size, int max_len) {
	int i, j, k, l, *len, *index, *bits, *link;

	len = (int *)alloc_mem(size * sizeof(int));
	index = (int *)alloc_mem(size * sizeof(int));
	bits = (int *)alloc_mem(size * sizeof(int));
	link = (int *)alloc_mem(size * sizeof(int));
	for (i = 0; i < size; i++) {
		len[i] = 0;
		index[i] = i;
		link[i] = -1;
	}
	/* sort in decreasing order of frequency */
	for (i = size -1; i > 0; i--) {
		for (j = 0; j < i; j++) {
			if (hist[index[j]] < hist[index[j + 1]]) {
				k = index[j + 1];
				index[j + 1] = index[j];
				index[j] = k;
			}
		}
	}
	for (i = 0; i < size; i++) {
		bits[i] = index[i];	/* reserv a sorted index table */
	}
	for (i = size - 1; i > 0; i--) {
		k = index[i - 1];
		l = index[i];
		hist[k] += hist[l];
		len[k]++;
		while (link[k] >= 0) {
			k = link[k];
			len[k]++;
		}
		link[k] = l;
		len[l]++;
		while (link[l] >= 0) {
			l = link[l];
			len[l]++;
		}
		for (j = i - 1; j > 0; j--) {
			if (hist[index[j - 1]] < hist[index[j]]) {
				k = index[j];
				index[j] = index[j - 1];
				index[j - 1] = k;
			} else {
				break;
			}
		}
	}
	/* limit the maximum code length to max_len */
	for (i = 0; i < size; i++) {
		index[i] = bits[i];	/* restore the index table */
		bits[i] = 0;
	}
	for (i = 0; i < size; i++) {
		bits[len[i]]++;
	}
	for (i = size - 1; i > max_len; i--) {
		while (bits[i] > 0) {
			j = i - 2;
			while(bits[j] == 0) j--;
			bits[i] -= 2;
			bits[i - 1]++;
			bits[j + 1] += 2;
			bits[j]--;
		}
	}
	for (i = k = 0; i < size; i++) {
		for (j = 0; j < bits[i]; j++) {
			len[index[k++]] = i;
		}
	}
	free(link);
	free(bits);
	free(index);
	return (len);
}

void gen_huffcode(VLC *vlc) {
	int i, j, *idx, *len;
	uint k;

	vlc->index = idx = (int *)alloc_mem(vlc->size * sizeof(int));
	vlc->off = (int *)alloc_mem(vlc->max_len * sizeof(int));
	vlc->code = (uint *)alloc_mem(vlc->size * sizeof(int));
	len = vlc->len;

	/* sort in increasing order of code length */
	for (i = 0; i < vlc->size; i++) {
		idx[i] = i;
	}

	for (i = vlc->size -1; i > 0; i--) {
		for (j = 0; j < i; j++) {
			if (len[idx[j]] > len[idx[j + 1]]) {
				k = idx[j + 1];
				idx[j + 1] = idx[j];
				idx[j] = k;
			}
		}
	}

	k = 0;

	for (j = 0; j < vlc->max_len; j++) {
		vlc->off[j] = -1;
	}

	j = len[idx[0]];

	for (i = 0; i < vlc->size; i++) {
		if (j < len[idx[i]]) {
			k <<= (len[idx[i]] - j);
			j = len[idx[i]];
		}

		vlc->code[idx[i]] = k++;
		vlc->off[j - 1] = i;
	}

	return;
}

VLC *make_vlc(uint *hist, int size, int max_len) {
	VLC *vlc;

	vlc = (VLC *)alloc_mem(sizeof(VLC));
	vlc->size = size;
	vlc->max_len = max_len;
	vlc->len = gen_hufflen(hist, size, max_len);
	gen_huffcode(vlc);
	return (vlc);
}

void free_vlc(VLC *vlc) {
	free(vlc->code);
	free(vlc->off);
	free(vlc->index);
	free(vlc->len);
	free(vlc);
	return;
}

VLC **init_vlcs(PMODEL ***pmodels, int num_group, int num_pmodel) {
	VLC **vlcs, *vlc;
	PMODEL *pm;
	int gr, k;

	vlcs = (VLC **)alloc_2d_array(num_group, num_pmodel, sizeof(VLC));
	for (gr = 0; gr < num_group; gr++) {
		for (k = 0; k < num_pmodel; k++) {
			vlc = &vlcs[gr][k];
			pm = pmodels[gr][k];
			vlc->size = pm->size;
			vlc->max_len = VLC_MAXLEN;
			vlc->len = gen_hufflen(pm->freq, pm->size, VLC_MAXLEN);
			gen_huffcode(vlc);
		}
	}
	return (vlcs);
}

/*
  Natural logarithm of the gamma function
  cf. "Numerical Recipes in C", 6.1
  http://www.ulib.org/webRoot/Books/Numerical_Recipes/bookcpdf.html
 */
double lngamma(double xx) {
	int j;
	double x,y,tmp,ser;
	double cof[6] = {
			76.18009172947146,	-86.50532032941677,
			24.01409824083091,	-1.231739572450155,
			0.1208650973866179e-2,	-0.5395239384953e-5
	};

	y = x = xx;
	tmp = x + 5.5 - (x + 0.5) * log(x + 5.5);
	ser = 1.000000000190015;
	for (j=0;j<=5;j++)
		ser += (cof[j] / ++y);
	return (log(2.5066282746310005 * ser / x) - tmp);
}

void set_freqtable(PMODEL *pm, double *pdfsamp, int num_subpm, int num_pmodel, int center, int idx, double sigma) {
	double shape, beta, norm, sw, x;
	int i, j, n;

	if (center == 0) sigma *= 2.0;

	if (idx < 0) {
		shape = 2.0;
	}
	else {
		shape = 3.2 * (idx + 1) / (double)num_pmodel;
	}

	/* Generalized Gaussian distribution */
	beta = exp(0.5*(lngamma(3.0/shape)-lngamma(1.0/shape))) / sigma;
	sw = 1.0 / (double)num_subpm;
	n = pm->size * num_subpm;
	center *= num_subpm;

	if (center == 0) {    /* one-sided distribution */
		for (i = 0; i < n; i++) {
			x = (double)i * sw;
			pdfsamp[i] = exp(-pow(beta * x, shape));
		}
	}
	else {
		for (i = center; i < n; i++) {
			x = (double)(i - (double)center + 0.5) * sw;
			pdfsamp[i + 1] = exp(-pow(beta * x, shape));
		}

		for (i = 0; i <= center; i++) {
			pdfsamp[center - i] =  pdfsamp[center + i + 1];
		}

		for (i = 0; i < n; i++) {
			if (i == center) {
				pdfsamp[i] = (2.0 + pdfsamp[i] + pdfsamp[i + 1]) / 2.0;
			}
			else {
				pdfsamp[i] = pdfsamp[i] + pdfsamp[i + 1];
			}
		}
	}

	for (j = 0; j < num_subpm; j++) {
		norm = 0.0;
		for (i = 0; i < pm->size; i++) {
			norm += pdfsamp[i * num_subpm + j];
		}

		norm = (double)(MAX_TOTFREQ - pm->size * MIN_FREQ) / norm;
		norm += 1E-8;	/* to avoid machine dependent rounding errors */
		pm->cumfreq[0] = 0;

		for (i = 0; i < pm->size; i++) {
			pm->freq[i] = norm * pdfsamp[i * num_subpm + j] + MIN_FREQ;
			pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
		}

		pm++;
	}
	return;
}

/*---------------------------- init_pmodels ------------------------*
 |  Function init_pmodels
 |
 |  Purpose:  Initializes the probability models
 |
 |  Parameters:
 |		num_group		--> Number of groups ???? (IN)
 |		num_pmodel		--> Number of probability models (IN)
 |		pm_accuracy		--> Probability model accuracy (IN)
 |		pm_idx			--> ?? --> NULL (IN)
 |		sigma			--> Predefined vectors (IN)
 |		size			--> Image maximum value + 1 (IN)
 |
 |  Returns:  PMODEL*** --> returns a PMODEL type structure
 *-------------------------------------------------------------------*/
PMODEL ***init_pmodels(int num_group, int num_pmodel, int pm_accuracy, int *pm_idx, double *sigma, int size) {
	PMODEL ***pmodels, *pmbuf, *pm;
	int gr, i, j, num_subpm, num_pm, ssize, idx;
	double *pdfsamp;

	//??
	if (pm_accuracy < 0) {
		num_subpm = 1;
		ssize = 1;
	}
	else {
		num_subpm = 1 << pm_accuracy;
		ssize = size;
		size = size + ssize - 1;
	}

	//Defines the number of probability models
	num_pm = (pm_idx != NULL)? 1 : num_pmodel;

	pmodels = (PMODEL ***)alloc_2d_array(num_group, num_pm, sizeof(PMODEL *));
	pmbuf = (PMODEL *)alloc_mem(num_group * num_pm * num_subpm * sizeof(PMODEL));

	//Vector initializations
	for (gr = 0; gr < num_group; gr++) {
		for (i = 0; i < num_pm; i++) {
			pmodels[gr][i] = pmbuf;

			for (j = 0; j < num_subpm; j++) {
				pm = pmbuf++;
				pm->id = i;
				pm->size = size;
				pm->freq = (uint *)alloc_mem((size * 2 + 1) * sizeof(uint));
				pm->cumfreq = &pm->freq[size];

				if (pm_idx == NULL) {
					pm->cost = (float *)alloc_mem((size + ssize) * sizeof(float));
					pm->subcost = &pm->cost[size];
				}
			}
		}
	}

	pdfsamp = alloc_mem((size * num_subpm + 1) * sizeof(double));

	for (gr = 0; gr < num_group; gr++) {
		for (i = 0; i < num_pm; i++) {
			if (pm_idx != NULL) {
				idx = pm_idx[gr];
			}
			else if (num_pm > 1) {
				idx = i;
			}
			else {
				idx = -1;
			}

			set_freqtable(pmodels[gr][i], pdfsamp, num_subpm, num_pmodel, ssize - 1, idx, sigma[gr]);
		}
	}

	free(pdfsamp);

	return (pmodels);
}

/* probability model for coefficients and thresholds */
void set_spmodel(PMODEL *pm, int size, int m) {
	int i, sum;
	double p;

	pm->size = size;

	if (m >= 0) {
		p = 1.0 / (double)(1 << (m % 8));
		sum = 0;

		for (i = 0; i < pm->size; i++) {
			pm->freq[i] = exp(-p * i) * (1 << 10);

			if (pm->freq[i] == 0) pm->freq[i]++;

			sum += pm->freq[i];
		}

		if (m & 8) pm->freq[0] = (sum - pm->freq[0]);	/* weight for zero */
	}
	else {
		for (i = 0; i < pm->size; i++) {
			pm->freq[i] = 1;
		}
	}

	pm->cumfreq[0] = 0;

	for (i = 0; i < pm->size; i++) {
		pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
	}

	return;
}

int *init_ctx_weight(int prd_order, int back_prd_order, int for_prd_order, int delta) {
	int *ctx_weight, k;
	double dy, dx;

	ctx_weight = (int *)alloc_mem((prd_order + back_prd_order + for_prd_order) * sizeof(int));

	for (k = 0; k < prd_order; k++) {
		dy = dyx[k].y;
		dx = dyx[k].x;

		ctx_weight[k] = 64.0 / sqrt(dy * dy + dx * dx) + 0.5;
	}

	if (back_prd_order > 0) {
		ctx_weight[prd_order] = 64.0 / sqrt(delta * delta) + 0.5;

		for (k = 0; k < back_prd_order - 1; k++) {
			dy = idyx[k].y;
			dx = idyx[k].x;

			ctx_weight[k + prd_order + 1] = 64.0 / sqrt(delta * delta + dy * dy + dx * dx) + 0.5;
		}
	}

	if (for_prd_order > 0) {
		ctx_weight[prd_order + back_prd_order] = 64.0 / sqrt(delta * delta) + 0.5;

		for (k = 0; k < for_prd_order - 1; k++) {
			dy = idyx[k].y;
			dx = idyx[k].x;

			ctx_weight[k + prd_order + back_prd_order + 1] = 64.0 / sqrt(delta * delta + dy * dy + dx * dx) + 0.5;
		}
	}

	return (ctx_weight);
}

int e2E(int e, int prd, int flag, int maxval) {
	int E, th;

	E = (e > 0)? e : -e;
	th = (prd < ((maxval + 1) >> 1))? prd : maxval - prd;

	if (E > th) {
		E += th;
	}
	else if (flag) {
		E = (e < 0)? (E << 1) - 1 : (E << 1);
	}
	else {
		E = (e > 0)? (E << 1) - 1 : (E << 1);
	}

	return (E);
}

int E2e(int E, int prd, int flag, int maxval) {
	int e, th;

	th = (prd < ((maxval + 1) >> 1))? prd : maxval - prd;

	if (E > (th << 1)) {
		e = (prd < ((maxval + 1) >> 1))? E - th : th - E;
	} else if (flag) {
		e = (E & 1)? -((E >> 1) + 1) : (E >> 1);
	} else {
		e = (E & 1)? (E >> 1) + 1 : -(E >> 1);
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
			ref[0] = ref[1] = ref[2] = class[y][x-1];
		}
	}
	else {
		ref[0] = class[y-1][x];
		ref[1] = (x == 0)? class[y-1][x] : class[y][x-1];
		ref[2] = (x + bsize >= width)? class[y-1][x] : class[y-1][x+bsize];
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

	return;
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
		dif = (unsigned)cur - prev;
	}
	prev = cur;

#ifdef HAVE_CLOCK
	return ((double)dif / CLOCKS_PER_SEC);
#else
	return ((double)dif / CLK_TCK);
#endif
}

int **select_bref(int bframes) {
	int (*aux)[5] = NULL;

	switch(bframes) {
		case 2:
			aux = bref1;
			break;
		case 3:
			aux = bref2;
			break;
		case 4:
			aux = bref3;
			break;
		case 5:
			aux = bref4;
			break;
		case 6:
			aux = bref5;
			break;
		case 7:
			aux = bref6;
			break;
		case 8:
			aux = bref7;
			break;
		case 9:
			aux = bref8;
			break;
		case 10:
			aux = bref9;
			break;
		default:
			aux = bref3;
			break;
	}

	return aux;
}
