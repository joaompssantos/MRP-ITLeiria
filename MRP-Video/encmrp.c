#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "mrp.h"

extern POINT dyx[];
extern POINT idyx[];
extern double sigma_h[], sigma_a[];
extern double qtree_prob[];

/*--------------------------- init_ref_offset -----------------------*
 |  Function init_ref_offset
 |
 |  Purpose:  Calculates the offset of a reference to a given pixel
 |
 |  Parameters:
 |		img				--> Image structure (IN)
 |		prd_order		--> Number of intra reference pixels (IN)
 |		inter_prd_order	--> Number of inter reference pixels (IN)
 |
 |  Returns:  int*** --> Array with the references offset
 *-------------------------------------------------------------------*/
int ***init_ref_offset(IMAGE *img, int prd_order, int back_prd_order, int for_prd_order) {
	int ***roff, *ptr;
	int x, y, dx, dy, k;
	int min_dx, max_dx, min_dy, max_dy;
	int bmin_dx, bmax_dx, bmin_dy, bmax_dy;
	int fmin_dx, fmax_dx, fmin_dy, fmax_dy;
	int min_abs_dx, max_abs_dx, min_abs_dy, max_abs_dy;
	int full_prd_order = prd_order + back_prd_order + for_prd_order;

	min_dx = max_dx = min_dy = max_dy = 0;
	bmin_dx = bmax_dx = bmin_dy = bmax_dy = 0;
	fmin_dx = fmax_dx = fmin_dy = fmax_dy = 0;

	//Values to check for special cases
	for (k = 0; k < prd_order; k++) {
		dy = dyx[k].y;
		dx = dyx[k].x;

		if (dy < min_dy) min_dy = dy;
		if (dx < min_dx) min_dx = dx;
		if (dy > max_dy) max_dy = dy;
		if (dx > max_dx) max_dx = dx;
	}

	for (k = 0; k < back_prd_order - 1; k++) {
		dy = idyx[k].y;
		dx = idyx[k].x;

		if (dy < bmin_dy) bmin_dy = dy;
		if (dy > bmax_dy) bmax_dy = dy;
		if (dx < bmin_dx) bmin_dx = dx;
		if (dx > bmax_dx) bmax_dx = dx;
	}

	for (k = 0; k < for_prd_order - 1; k++) {
		dy = idyx[k].y;
		dx = idyx[k].x;

		if (dy < fmin_dy) fmin_dy = dy;
		if (dy > fmax_dy) fmax_dy = dy;
		if (dx < fmin_dx) fmin_dx = dx;
		if (dx > fmax_dx) fmax_dx = dx;
	}

	min_abs_dy = (min_dy < bmin_dy && min_dy < fmin_dy ? min_dy : bmin_dy < fmin_dy ? bmin_dy : fmin_dy);
	min_abs_dx = (min_dx < bmin_dx && min_dx < fmin_dx ? min_dx : bmin_dx < fmin_dx ? bmin_dx : fmin_dx);

	if (back_prd_order <= 1 && for_prd_order <= 1) {
		max_abs_dy = 0;
	}
	else {
		max_abs_dy = (fmax_dy > bmax_dy ? fmax_dy : bmax_dy);
	}

	max_abs_dx = (max_dx > bmax_dx && max_dx > fmax_dx ? max_dx : bmax_dx > fmax_dx ? bmax_dx : fmax_dx);

	roff = (int ***)alloc_2d_array(img->height, img->width, sizeof(int *));

	//Cycle that runs for all the pixels
	for (y = 0; y < img->height; y++) {
		for (x = 0; x < img->width; x++) {
			ptr = (int *) alloc_mem((full_prd_order) * sizeof(int));
			//Conditions to check which references are available for each pixel
			if (y == 0) {
				if (x == 0) {
					roff[y][x] = ptr;
					dx = 0;
					dy = img->height;

					for (k = 0; k < prd_order; k++) {
						*ptr++ = dy * img->width + dx; //Points to a line filled with 128 (if max_val = 256)
					}
				}
				else if (x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
					roff[y][x] = ptr;
					dy = 0;

					for (k = 0; k < prd_order; k++) {
						dx = dyx[k].x;

						if (x + dx < 0) dx = -x;
						else if (dx >= 0) dx = -1;

						*ptr++ = dy * img->width + dx;
					}
				}
				else {
					roff[y][x] = roff[y][x - 1];
					free(ptr);
				}
			}
			else if (y + min_abs_dy <= 0) {
				if (x == 0) {
					roff[y][x] = ptr;

					for (k = 0; k < prd_order; k++) {
						dy = dyx[k].y;

						if (y + dy < 0) dy = -y;
						else if (dy >= 0) dy = -1;

						dx = dyx[k].x;

						if (x + dx < 0) dx = -x;

						*ptr++ = dy * img->width + dx;
					}
				}
				else if (x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
					roff[y][x] = ptr;

					for (k = 0; k < prd_order; k++) {
						dy = dyx[k].y;

						if (y + dy < 0) dy = -y;
						dx = dyx[k].x;

						if (x + dx < 0) dx = -x;
						else if (x + dx >= img->width) {
							dx = img->width - x - 1;
						}

						*ptr++ = dy * img->width + dx;
					}
				}
				else {
					roff[y][x] = roff[y][x - 1];
					free(ptr);
				}
			}
			else if (y + max_abs_dy >= img->height) {
				if (x == 0 || x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
					roff[y][x] = ptr;

					for (k = 0; k < prd_order; k++) {
						*ptr++ = roff[y - 1][x][k];
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

			//Back reference offset
			if (back_prd_order != 0) {
				int base = -img->width * (img->height + 1);

				if (y == 0) {
					if (x == 0) {
						*ptr++ = base;

						for (k = 0; k < back_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if (y + dy < 0 || x + dx < 0) {
								*ptr++ = base;
							}
							else {
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else if (x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
						*ptr++ = base;

						for (k = 0; k < back_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if (y + dy < 0 || x + dx < 0 || x + dx >= img->width) {
								*ptr++ = base;
							}
							else {
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else {
						roff[y][x] = roff[y][x - 1];
					}
				}
				else if (y + min_abs_dy <= 0) {
					if (x == 0) {
						*ptr++ = base;

						for (k = 0; k < back_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if (y + dy < 0 || x + dx < 0) {
								*ptr++ = base;
							}
							else {
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else if (x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
						*ptr++ = base;

						for (k = 0; k < back_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if (y + dy < 0 || x + dx < 0 || x + dx >= img->width) {
								*ptr++ = base;
							}
							else {
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else {
						roff[y][x] = roff[y][x - 1];
					}
				}
				else if (y + max_abs_dy >= img->height) {
					if (x == 0 || x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
						*ptr++ = base;

						for (k = 0; k < back_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if (y + dy >= img->height || x + dx < 0 || x + dx >= img->width) {
								*ptr++ = base;
							}
							else {
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else {
						roff[y][x] = roff[y][x - 1];
					}
				}
				else {
					roff[y][x] = roff[y - 1][x];
				}
			}

			//Forward reference offset
			if (for_prd_order != 0) {
				int base = img->width * (img->height + 1);

				if (y == 0) {
					if (x == 0) {
						*ptr++ = base;

						for (k = 0; k < for_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if (y + dy < 0 || x + dx < 0) {
								*ptr++ = base;
							}
							else {
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else if (x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
						*ptr++ = base;

						for (k = 0; k < for_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if (y + dy < 0 || x + dx < 0 || x + dx >= img->width) {
								*ptr++ = base;
							}
							else {
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else {
						roff[y][x] = roff[y][x - 1];
					}
				}
				else if (y + min_abs_dy <= 0) {
					if (x == 0) {
						*ptr++ = base;

						for (k = 0; k < for_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if (y + dy < 0 || x + dx < 0) {
								*ptr++ = base;
							}
							else {
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else if (x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
						*ptr++ = base;

						for (k = 0; k < for_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if (y + dy < 0 || x + dx < 0 || x + dx >= img->width) {
								*ptr++ = base;
							}
							else {
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else {
						roff[y][x] = roff[y][x - 1];
					}
				}
				else if (y + max_abs_dy >= img->height) {
					if (x == 0 || x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
						*ptr++ = base;

						for (k = 0; k < for_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if (y + dy >= img->height || x + dx < 0 || x + dx >= img->width) {
								*ptr++ = base;
							}
							else {
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else {
						roff[y][x] = roff[y][x - 1];
					}
				}
				else {
					roff[y][x] = roff[y - 1][x];
				}
			}
		}
	}

	return (roff);
}

/*---------------------------- init_encoder ------------------------*
 |  Function init_encoder
 |
 |  Purpose:  Initializes the encoder structure
 |
 |  Parameters:
 |      img		 		--> Video/Image file structure (IN)
 |		num_class		--> Number of classes to use (number of different predictors) (IN)
 |		num_group		--> Number of groups ???? (IN)
 |		prd_order		--> Order of the predictors (number of pixels to use) (IN)
 |		inter_prd_order	--> Order of the inter predictors (number of pixels to use) (IN)
 |		coef_precision	--> Precision of the coefficients (IN)
 |		f_huffman		--> Huffman coding flag (IN)
 |		quadtree_depth	--> Quadtree flag (IN)
 |		num_pmodel		--> Number of probability models (IN)
 |		pm_accuracy		--> Probability model accuracy (IN)
 |
 |  Returns:  ENCODER* --> returns a encoder type structure
 *-------------------------------------------------------------------*/
ENCODER *init_encoder(IMAGE *img, IMAGE *back_ref_img, IMAGE *for_ref_img, int **back_ref_error, int **for_ref_error, int num_class, int num_group, int prd_order, int back_prd_order, int for_prd_order, int coef_precision, int f_huffman, int quadtree_depth, int num_pmodel, int pm_accuracy, int delta) {
	//Declare new encoder struct
	ENCODER *enc;
	int x, y, i, j;
	double c;
	int full_prd_order;

	//Allocation of the memory for the encoder struct
	enc = (ENCODER *)alloc_mem(sizeof(ENCODER));

	//Copy of the video/image properties to encoder
	enc->height = img->height;
	enc->width = img->width;
	enc->maxval = img->maxval;

	// Copy of the encoding parameters
	enc->num_class = num_class;
	enc->num_group = num_group; // ?? (= 16)
	enc->prd_order = prd_order; // K
	enc->back_prd_order = back_prd_order;
	enc->for_prd_order = for_prd_order;
	enc->coef_precision = coef_precision; // P
	enc->f_huffman = f_huffman; // -h
	enc->quadtree_depth = quadtree_depth; // -f
	enc->num_pmodel = num_pmodel; // V
	enc->pm_accuracy = pm_accuracy; // A
	enc->maxprd = enc->maxval << enc->coef_precision; // ??
	enc->delta = delta;
	full_prd_order = enc->prd_order + enc->back_prd_order + enc->for_prd_order;

	// Alloc memory to predictors array
	enc->predictor = (int **)alloc_2d_array(enc->num_class, full_prd_order, sizeof(int));

	// Alloc memory to ??
	enc->th = (int **)alloc_2d_array(enc->num_class, enc->num_group, sizeof(int));

	// enc->th initializing cycle
	for (i = 0; i < enc->num_class; i++) {
		for (j = 0; j < enc->num_group - 1; j++) {
			enc->th[i][j] = 0;
		}
		enc->th[i][enc->num_group - 1] = MAX_UPARA + 1;
	}

	// More memory allocation
	enc->upara = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));
	enc->prd = (int **)alloc_2d_array(enc->height, enc->width, sizeof(int));
	// Memory allocation for the original image
	enc->org = (int ***)alloc_3d_array(enc->height + 1, enc->width, 3, sizeof(int));
	// Memory allocation for the error image ??
	enc->err = (int ***)alloc_3d_array(enc->height + 1, enc->width, 3, sizeof(int));

	// Quadtree map what?
	if (enc->quadtree_depth > 0) {
		y = (enc->height + MAX_BSIZE - 1) / MAX_BSIZE;
		x = (enc->width + MAX_BSIZE - 1) / MAX_BSIZE;
		for (i = enc->quadtree_depth - 1; i >= 0; i--) {
			enc->qtmap[i] = (char **)alloc_2d_array(y, x, sizeof(char));
			y <<= 1;
			x <<= 1;
		}
	}

	enc->ctx_weight = init_ctx_weight(enc->prd_order, enc->back_prd_order, enc->for_prd_order, enc->delta);

	// Class and group arrays
	enc->class = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));
	enc->group = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));

	// Copy original images to encoder struct, initialize group array
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			enc->group[y][x] = 0;
			enc->org[1][y][x] = img->val[y][x];
		}
	}

	if (back_ref_img != NULL) {
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				enc->org[0][y][x] = back_ref_img->val[y][x];
			}
		}
	}

	if (for_ref_img != NULL) {
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				enc->org[2][y][x] = for_ref_img->val[y][x];
			}
		}
	}

	if (back_ref_error != NULL) {
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				enc->err[0][y][x] = back_ref_error[y][x];
			}
		}
	}

	if (for_ref_error != NULL) {
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				enc->err[2][y][x] = for_ref_error[y][x];
			}
		}
	}

	// Auxiliary values
	enc->org[2][enc->height][0] = (enc->maxval + 1) >> 1;
	enc->org[1][enc->height][0] = (enc->maxval + 1) >> 1;
	enc->org[0][enc->height][0] = (enc->maxval + 1) >> 1;
	enc->err[2][enc->height][0] = (enc->maxval + 1) >> 2;
	enc->err[1][enc->height][0] = (enc->maxval + 1) >> 2;
	enc->err[0][enc->height][0] = (enc->maxval + 1) >> 2;

	// Table used for the quantization of the context variable
	enc->uquant = (char **)alloc_2d_array(enc->num_class, MAX_UPARA + 1, sizeof(char));

	for (i = 0; i < enc->num_class; i++) {
		for (j = 0; j <= MAX_UPARA; j++) {
			enc->uquant[i][j] = enc->num_group - 1;
		}
	}

	enc->econv = (int **)alloc_2d_array(enc->maxval+1, (enc->maxval<<1) + 1, sizeof(int));
	enc->bconv = (img_t *)alloc_mem((enc->maxprd + 1) * sizeof(img_t));
	enc->fconv = (img_t *)alloc_mem((enc->maxprd + 1) * sizeof(img_t));

	enc->pmlist = (PMODEL **)alloc_mem(enc->num_group * sizeof(PMODEL *));

	enc->spm.freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
	enc->spm.cumfreq = &(enc->spm.freq[MAX_SYMBOL]);

	// Huffman coding
	if (enc->f_huffman == 1) {
		enc->sigma = sigma_h;
	} else {
		enc->sigma = sigma_a;
	}

	// ?
	enc->mtfbuf = (int *)alloc_mem(enc->num_class * sizeof(int));
	enc->coef_m = (int *)alloc_mem((full_prd_order) * sizeof(int));

	for (i = 0; i < full_prd_order; i++) {
		enc->coef_m[i] = 0;
	}

	enc->coef_cost = (cost_t **)alloc_2d_array(16, MAX_COEF + 1, sizeof(cost_t));

	for (i = 0; i < 16; i++) {
#ifdef OPT_SIDEINFO
		if (enc->f_huffman == 1) {
			for (j = 0; j <= MAX_COEF; j++) {
				enc->coef_cost[i][j] = ((j >> i) + i + 1);
				if (j > 0) enc->coef_cost[i][j] += 1.0;
			}
		} else {
			double p;
			set_spmodel(&enc->spm, MAX_COEF + 1, i);
			p = log(enc->spm.cumfreq[MAX_COEF + 1]);
			for (j = 0; j <= MAX_COEF; j++) {
				enc->coef_cost[i][j] = (p - log(enc->spm.freq[j])) / log(2.0);
				if (j > 0) enc->coef_cost[i][j] += 1.0;
			}
		}
#else
		for (j = 0; j <= MAX_COEF; j++) {
			enc->coef_cost[i][j] = 0;
		}
#endif
	}

	enc->th_cost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
	for (i = 0; i < MAX_UPARA + 2; i++) {
		enc->th_cost[i] = 0;
	}

	enc->class_cost = (cost_t *)alloc_mem(enc->num_class * sizeof(cost_t));

	c = log((double)enc->num_class) / log(2.0);

	for (i = 0; i < enc->num_class; i++) {
		enc->class_cost[i] = c;
	}
	for (i = 0; i < (QUADTREE_DEPTH << 3); i++) {
		enc->qtflag_cost[i] = 1.0;
	}
	return (enc);
}

int **get_enc_err(ENCODER *enc, int pos) {
	int y, x;
	int **error = (int **)alloc_2d_array(enc->height + 1, enc->width, sizeof(int));

	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			error[y][x] = enc->err[pos][y][x];
		}
	}

	return error;
}

void free_encoder(ENCODER *enc) {
	int dy, dx, x, y, i, j, k;
	int min_dx, max_dx, min_dy, max_dy;
	int bmin_dx, bmax_dx, bmin_dy, bmax_dy;
	int fmin_dx, fmax_dx, fmin_dy, fmax_dy;
	int min_abs_dx, max_abs_dx, min_abs_dy, max_abs_dy;
	int gr;
	int num_subpm;

	free(enc->predictor);
	free(enc->th);
	free(enc->upara);
	free(enc->prd);
	free(enc->org);
	free(enc->err);

	if (enc->quadtree_depth > 0) {
		for (x = enc->quadtree_depth - 1; x >= 0; x--) {
			free(enc->qtmap[x]);
		}
	}

	free(enc->class);
	free(enc->group);
	free(enc->uquant);
	free(enc->econv);
	free(enc->bconv);
	free(enc->fconv);
	free(enc->pmlist);
	free(enc->spm.freq);
	free(enc->mtfbuf);
	free(enc->coef_m);
	free(enc->coef_cost);
	free(enc->th_cost);
	free(enc->class_cost);
	free(enc->ctx_weight);

	if (enc->f_huffman == 0) {
		free(enc->rc);
	}

	min_dx = max_dx = min_dy = max_dy = 0;
	bmin_dx = bmax_dx = bmin_dy = bmax_dy = 0;
	fmin_dx = fmax_dx = fmin_dy = fmax_dy = 0;

	//Values to check for special cases
	for (k = 0; k < enc->prd_order; k++) {
		dy = dyx[k].y;
		dx = dyx[k].x;

		if (dy < min_dy) min_dy = dy;
		if (dx < min_dx) min_dx = dx;
		if (dy > max_dy) max_dy = dy;
		if (dx > max_dx) max_dx = dx;
	}

	for (k = 0; k < enc->back_prd_order - 1; k++) {
		dy = idyx[k].y;
		dx = idyx[k].x;

		if (dy < bmin_dy) bmin_dy = dy;
		if (dy > bmax_dy) bmax_dy = dy;
		if (dx < bmin_dx) bmin_dx = dx;
		if (dx > bmax_dx) bmax_dx = dx;
	}

	for (k = 0; k < enc->for_prd_order - 1; k++) {
		dy = idyx[k].y;
		dx = idyx[k].x;

		if (dy < fmin_dy) fmin_dy = dy;
		if (dy > fmax_dy) fmax_dy = dy;
		if (dx < fmin_dx) fmin_dx = dx;
		if (dx > fmax_dx) fmax_dx = dx;
	}

	min_abs_dy = (min_dy < bmin_dy && min_dy < fmin_dy ? min_dy : bmin_dy < fmin_dy ? bmin_dy : fmin_dy);
	min_abs_dx = (min_dx < bmin_dx && min_dx < fmin_dx ? min_dx : bmin_dx < fmin_dx ? bmin_dx : fmin_dx);

	if (enc->back_prd_order <= 1 && enc->for_prd_order <= 1) {
		max_abs_dy = 0;
	}
	else {
		max_abs_dy = (fmax_dy > bmax_dy ? fmax_dy : bmax_dy);
	}

	max_abs_dx = (max_dx > bmax_dx && max_dx > fmax_dx ? max_dx : bmax_dx > fmax_dx ? bmax_dx : fmax_dx);

	//Cycle that runs for all the pixels
	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			//Conditions to check which references are available for each pixel
			if ((y == 0 && (x == 0 || x + min_abs_dx <= 0 || x + max_abs_dx >= enc->width)) || \
				(y + min_abs_dy <= 0 && (x == 0 || x + min_abs_dx <= 0 || x + max_abs_dx >= enc->width)) || \
				(y + max_abs_dy >= enc->height && (x == 0 || x + min_abs_dx <= 0 || x + max_abs_dx >= enc->width))) {
				free(enc->roff[y][x]);
			}
		}
	}
	safefree((void **)&(enc->roff));

	if (enc->pm_accuracy < 0) {
		num_subpm = 1;
	}
	else {
		num_subpm = 1 << enc->pm_accuracy;
	}
	for (gr = 0; gr < enc->num_group; gr++) {
		for (i = 0; i < enc->num_pmodel; i++) {
			for (j = 0; j < num_subpm; j++) {
				PMODEL *aux = &enc->pmodels[gr][i][j];
				free(aux->freq);
				free(aux->cost);
			}
		}
	}
	free(enc->pmodels[0][0]);
	free(enc->pmodels);

	if (enc->f_huffman == 1) {
		int gr;
		VLC *aux;
		for (gr = 0; gr < enc->num_group; gr++) {
			for (k = 0; k < enc->num_pmodel; k++) {
				aux = &enc->vlcs[gr][k];
				free(aux->code);
				free(aux->off);
				free(aux->index);
				free(aux->len);
			}
		}
		free(enc->vlcs);
	}

	free(enc);
}

/*------------------------------ init_class -------------------------*
 |  Function init_class
 |
 |  Purpose:  Sets the class number of each pixel in a given frame
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |		frame			--> Frame to be processed (IN)
 |
 |  Returns:  void
 *-------------------------------------------------------------------*/
void init_class(ENCODER *enc) {
	int k, x, y, i, j, v, cl, sum, num_block;
	int *var, *tmp, **ptr;

	// Number of blocks in the frame
	num_block = enc->height * enc->width / (BASE_BSIZE * BASE_BSIZE);

	var = (int *)alloc_mem(num_block * sizeof(int));
	ptr = (int **)alloc_mem(num_block * sizeof(int *));

	//Variance calculation for each block
	for (k = 0; k < num_block; k++) {
		//Gives the position of the top right pixel of each block
		y = (k / (enc->width / BASE_BSIZE)) * BASE_BSIZE;
		x = (k % (enc->width / BASE_BSIZE)) * BASE_BSIZE;

		var[k] = sum = 0;

		//Run each pixel in a block
		for (i = 0; i < BASE_BSIZE; i++) {
			for (j = 0; j < BASE_BSIZE; j++) {
				v = enc->org[1][y + i][x + j];
				sum += v;
				var[k] += v * v;
			}
		}

		//Final result of the variance for one block
		var[k] -= sum * sum / (BASE_BSIZE * BASE_BSIZE);
		ptr[k] = &(var[k]);
	}

	//Variance sorting, from lowest to highest
	for (i = num_block - 1; i > 0; i--) {
		for (j = 0; j < i; j++) {
			if (*ptr[j] > *ptr[j + 1]) {
				tmp = ptr[j];
				ptr[j] = ptr[j + 1];
				ptr[j + 1] = tmp;
			}
		}
	}

	for (k = 0; k < num_block; k++) {
		//Calculates the number of the class for a block, given its position (1, 2, ...)
		cl = (k * enc->num_class) / num_block;
		//Determines the new position of a block after being sorted
		v = (int)(ptr[k] - var);

		//Gives the position of the top right pixel of each of the sorted blocks
		y = (v / (enc->width / BASE_BSIZE)) * BASE_BSIZE;
		x = (v % (enc->width / BASE_BSIZE)) * BASE_BSIZE;

		//Sets the class number for each pixel
		for (i = 0; i < BASE_BSIZE; i++) {
			for (j = 0; j < BASE_BSIZE; j++) {
				enc->class[y + i][x + j] = cl;
			}
		}
	}

	free(ptr);
	free(var);
}

void set_cost_model(ENCODER *enc, int f_mmse) {
	int gr, i, j, k;
	double a, b, c, var;
	PMODEL *pm;

	for (i = 0; i <= enc->maxval; i++) {
		for (j = 0; j <= (enc->maxval << 1); j++) {
			k = (i << 1) - j;
			enc->econv[i][j] = (k > 0)? (k - 1) : (-k);
		}
	}

	enc->encval = enc->err[1];

	for (gr = 0; gr < enc->num_group; gr++) {
		var = enc->sigma[gr] * enc->sigma[gr];

		if (f_mmse) {
			a = 0;
			b = 1.0;
		}
		else {
			a = 0.5 * log(2 * M_PI * var) / log(2.0);
			b = 1.0 / (2.0 * log(2.0) * var);
		}

		enc->pmlist[gr] = pm = enc->pmodels[gr][enc->num_pmodel >> 1];

		for (k = 0; k <= pm->size; k++) {
			c = (double)k * 0.5 + 0.25;
			pm->cost[k] = a + b * c * c;
		}

		pm->subcost[0] = 0.0;
	}

	for (k = 0; k <= enc->maxprd; k++) {
		enc->bconv[k] = 0;
		enc->fconv[k] = 0;
	}

	return;
}

void set_cost_rate(ENCODER *enc) {
	int gr, k, i, j, mask, shift, num_spm;
	double a, c;
	PMODEL *pm;

	if (enc->pm_accuracy < 0) {
		for (i = 0; i <= enc->maxval; i++) {
			for (j = 0; j <= (enc->maxval << 1); j++) {
				k = (j + 1) >> 1;
				enc->econv[i][j] = e2E(i - k, k, j&1, enc->maxval);
			}
		}
	}

	if (enc->pm_accuracy < 0) {
		num_spm = 1;
	}
	else {
		enc->encval = enc->org[1];
		mask = (1 << enc->pm_accuracy) - 1;
		shift = enc->coef_precision - enc->pm_accuracy;

		for (k = 0; k <= enc->maxprd; k++) {
			i = (enc->maxprd - k + (1 << shift) / 2) >> shift;
			enc->fconv[k] = (i & mask);
			enc->bconv[k] = (i >> enc->pm_accuracy);
		}

		num_spm = 1 << enc->pm_accuracy;
	}

	a = 1.0 / log(2.0);

	for (gr = 0; gr < enc->num_group; gr++) {
		for (i = 0; i < enc->num_pmodel; i++) {
			pm = enc->pmodels[gr][i];

			if (enc->f_huffman == 1) {
				for (k = 0; k < pm->size; k++) {
					pm->cost[k] = enc->vlcs[gr][pm->id].len[k];
				}

				pm->subcost[0] = 0.0;
			}
			else if (enc->pm_accuracy < 0) {
				for (k = 0; k < pm->size; k++) {
					pm->cost[k] = -a * log(pm->freq[k]);
				}

				c = pm->cumfreq[enc->maxval + 1];
				pm->subcost[0] = a * log(c);
			}
			else {
				for (j = 0; j < num_spm; j++) {
					for (k = 0; k < pm->size; k++) {
						pm->cost[k] = -a * log(pm->freq[k]);
					}
					for (k = 0; k <= enc->maxval; k++) {
						c = pm->cumfreq[k + enc->maxval + 1] - pm->cumfreq[k];
						pm->subcost[k] = a * log(c);
					}

					pm++;
				}
			}
		}
	}
}

/*-------------------------- prediction_region ----------------------*
 |  Function prediction_region
 |
 |  Purpose: "Calculates the prediction of each pixel"
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |		frame			--> Frame to be processed (IN)
 |		tly				--> Top left Y boundary
 |		tlx				--> Top left X boundary
 |		bry				--> Bottom right Y boundary
 |		brx				--> Bottom right X boundary
 |
 |  Returns:  void
 *-------------------------------------------------------------------*/
void predict_region(ENCODER *enc, int tly, int tlx, int bry, int brx) {
	int x, y, k, cl, prd, org;
	int *coef_p;
	int *prd_p;
	int *roff_p, **roff_pp;
	int *err_p, *org_p;
	char *class_p;

	int full_prd_order = enc->prd_order + enc->back_prd_order + enc->for_prd_order;

	//Runs all the pixels in a frame (due to the used boundaries)
	for (y = tly; y < bry; y++) {
		class_p = &enc->class[y][tlx];
		org_p = &enc->org[1][y][tlx];
		roff_pp = &enc->roff[y][tlx];
		err_p = &enc->err[1][y][tlx];
		prd_p = &enc->prd[y][tlx];

		for (x = tlx; x < brx; x++) {
			cl = *class_p++;
			roff_p = *roff_pp++;
			coef_p = enc->predictor[cl];
			prd = 0;

			for (k = 0; k < full_prd_order; k++) {
				prd += org_p[*roff_p++] * (*coef_p++);
			}

			org = *org_p++;
			*prd_p++ = prd;

			//Boundaries of the predicted values
			if (prd < 0) prd = 0;
			else if (prd > enc->maxprd) prd = enc->maxprd;

			prd >>= (enc->coef_precision - 1);
			*err_p++ = enc->econv[org][prd];
		}
	}
}

int calc_uenc(ENCODER *enc, int y, int x) {
	int u, k, *err_p, *roff_p, *wt_p;
	int full_prd_order = enc->prd_order + enc->back_prd_order + enc->for_prd_order;

	err_p = &enc->err[1][y][x];
	roff_p = enc->roff[y][x];
	wt_p = enc->ctx_weight;
	u = 0;

	for (k = 0; k < full_prd_order; k++) {
		u += err_p[*roff_p++] * (*wt_p++);
	}

	u >>= 6;

	if (u > MAX_UPARA) u = MAX_UPARA;

	return (u);
}

cost_t calc_cost(ENCODER *enc, int tly, int tlx, int bry, int brx) {
	cost_t cost;
	int x, y, u, cl, gr, prd, e, base, frac;
	int *upara_p, *prd_p, *encval_p;
	char *class_p, *group_p;
	PMODEL *pm;

	//    bry += (UPEL_DIST - 1);
	//    tlx -= (UPEL_DIST - 1);
	//    brx += (UPEL_DIST - 1);
	if (bry > enc->height) bry = enc->height;
	if (tlx < 0) tlx = 0;
	if (brx > enc->width) brx = enc->width;

	cost = 0;

	for (y = tly; y < bry; y++) {
		class_p  = &enc->class[y][tlx];
		group_p  = &enc->group[y][tlx];
		upara_p  = &enc->upara[y][tlx];
		encval_p = &enc->encval[y][tlx];
		prd_p    = &enc->prd[y][tlx];

		for (x = tlx; x < brx; x++) {
			cl = *class_p++;
			*upara_p++ = u = calc_uenc(enc, y, x);
			*group_p++ = gr = enc->uquant[cl][u];
			e = *encval_p++;
			prd = *prd_p++;

			if (prd < 0) prd = 0;
			else if (prd > enc->maxprd) prd = enc->maxprd;

			base = enc->bconv[prd];
			frac = enc->fconv[prd];
			pm = enc->pmlist[gr] + frac;
			cost += pm->cost[base + e] + pm->subcost[base];
		}
	}

	return (cost);
}

/*-------------------------- design_predictor -----------------------*
 |  Function design_predictor
 |
 |  Purpose:  Calculates the predictor coefficients for each class
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |		frame			--> Frame to be processed (IN)
 |		f_mmse			--> Coding type (Huffman/Arithmetic) (IN)
 |
 |  Returns:  cost_t --> cost of the designed predictor
 *-------------------------------------------------------------------*/
cost_t design_predictor(ENCODER *enc, int f_mmse) {
	double **mat, *weight, w, e, d, pivot;
	int x, y, i, j, k, cl, gr, pivpos, *index, *roff_p, *org_p;

	int full_prd_order = enc->prd_order + enc->back_prd_order + enc->for_prd_order;

	mat = (double **)alloc_2d_array(full_prd_order, full_prd_order + 1, sizeof(double));
	index = (int *)alloc_mem(sizeof(int) * full_prd_order);
	weight = (double *)alloc_mem(sizeof(double) * enc->num_group);

	//Weight choice
	for (gr = 0; gr < enc->num_group; gr++) {
		if (f_mmse) {
			weight[gr] = 1.0;
		} else {
			weight[gr] = 1.0 / (enc->sigma[gr] * enc->sigma[gr]);
		}
	}

	//Cycle that runs each class in a given frame
	for (cl = 0; cl < enc->num_class; cl++) {
		//Variable mat initialization
		for (i = 0; i < full_prd_order; i++) {
			for (j = 0; j <= full_prd_order; j++) {
				mat[i][j] = 0.0;
			}
		}

		//Run each pixel
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				if (enc->class[y][x] != cl) { //Check if the pixel is member of the current class
					x += BASE_BSIZE - 1;
					continue;
				}

				gr = enc->group[y][x];
				roff_p = enc->roff[y][x]; //This variable has the position of the reference pixels
				org_p = &enc->org[1][y][x];

				//Fills the matrix mat with the reference pixels values
				for (i = 0; i < full_prd_order; i++) {
					w = weight[gr] * org_p[roff_p[i]];

					for (j = i; j < full_prd_order; j++) {
						mat[i][j] += w * org_p[roff_p[j]];
					}

					mat[i][full_prd_order] += w * org_p[0]; //Results column with the value of the pixel to be predicted
				}
			}
		}

		//Makes the bottom diagonal equal to the top diagonal
		for (i = 0; i < full_prd_order; i++) {
			index[i] = i;

			for (j = 0; j < i; j++) {
				mat[i][j] = mat[j][i];
			}
		}

		//"Gaussian elimination"
		for (i = 0; i < full_prd_order; i++) {
			pivpos = i;
			pivot = fabs(mat[index[i]][i]);

			//Sorts the matrix pivots
			for (k = i + 1; k < full_prd_order; k++) {
				if (fabs(mat[index[k]][i]) > pivot) {
					pivot = fabs(mat[index[k]][i]);
					pivpos = k;
				}
			}

			//Change the indexes
			k = index[i];
			index[i] = index[pivpos];
			index[pivpos] = k;

			//If the pivot is not zero the actual Gaussian elimination is performed
			if (pivot > 1E-10) {
				d = mat[index[i]][i];

				for (j = i; j <= full_prd_order; j++) {
					mat[index[i]][j] /= d;
				}

				for (k = 0; k < full_prd_order; k++) {
					if (k == i) continue;

					d = mat[index[k]][i];

					for (j = i; j <= full_prd_order; j++) {
						mat[index[k]][j] -= d * mat[index[i]][j];
					}
				}
			}
		}

		w = (1 << enc->coef_precision);
		e = 0.0;

		//Rounds the coefficients and stores the error
		for (i = 0; i < full_prd_order; i++) {
			if (fabs(mat[index[i]][i]) > 1E-10) { //Checks if a line is not zero
				d = mat[index[i]][full_prd_order] * w;
			}
			else {
				d = 0.0;
			}

			k = d;
			if (k > d) k--;

			if (k < -MAX_COEF) {
				k = d = -MAX_COEF;
			}
			else if (k > MAX_COEF) {
				k = d = MAX_COEF;
			}

			enc->predictor[cl][i] = k; //Coefficient for a given reference in the predictor
			d -= k;
			e += d;
			mat[index[i]][full_prd_order] = d; //Stores the error
		}

		/* minimize mean rounding errors */
		k = e + 0.5;
		for (;k > 0; k--) {
			d = 0;
			for (j = i = 0; i < full_prd_order; i++) {
				if (mat[index[i]][full_prd_order] > d) {
					d = mat[index[i]][full_prd_order];
					j = i;
				}
			}

			if (enc->predictor[cl][j] < MAX_COEF) enc->predictor[cl][j]++;

			mat[index[j]][full_prd_order] = 0;
		}
	}

	free(weight);
	free(index);
	free(mat);

	predict_region(enc, 0, 0, enc->height, enc->width);

	return (calc_cost(enc, 0, 0, enc->height, enc->width));
}

cost_t optimize_group(ENCODER *enc) {
	cost_t cost, min_cost, **cbuf, *dpcost, *cbuf_p, *thc_p;
	int x, y, th1, th0, k, u, cl, gr, prd, e, base, frac;
	int **trellis, *tre_p;
	PMODEL *pm, **pm_p;

	trellis = (int **)alloc_2d_array(enc->num_group, MAX_UPARA + 2, sizeof(int));
	dpcost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
	cbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2, sizeof(cost_t));
	thc_p = enc->th_cost;

	for (k = 0; k < MAX_UPARA + 2; k++) trellis[0][k] = 0;
	/* Dynamic programming */
	for (cl = 0; cl < enc->num_class; cl++) {
		for (gr = 0; gr < enc->num_group; gr++) {
			cbuf_p = cbuf[gr];

			for (u = 0; u < MAX_UPARA + 2; u++) {
				cbuf_p[u] = 0;
			}
		}

		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				if (enc->class[y][x] == cl) {
					u = enc->upara[y][x] + 1;
					e = enc->encval[y][x];
					prd = enc->prd[y][x];

					if (prd < 0) prd = 0;
					else if (prd > enc->maxprd) prd = enc->maxprd;

					base = enc->bconv[prd];
					frac = enc->fconv[prd];
					pm_p = enc->pmlist;

					for (gr = 0; gr < enc->num_group; gr++) {
						pm = (*pm_p++) + frac;
						cbuf[gr][u] += pm->cost[base + e] + pm->subcost[base];
					}
				}
			}
		}

		for (gr = 0; gr < enc->num_group; gr++) {
			cbuf_p = cbuf[gr];

			for (u = 1; u < MAX_UPARA + 2; u++) {
				cbuf_p[u] += cbuf_p[u - 1];
			}
		}

		cbuf_p = cbuf[0];

		for (u = 0; u < MAX_UPARA + 2; u++) {
			dpcost[u] = cbuf_p[u];
		}

		for (gr = 1; gr < enc->num_group; gr++) {
			cbuf_p = cbuf[gr];
			tre_p = trellis[gr - 1];

			/* minimize (cbuf_p[th1] - cbuf_p[th0] + dpcost[th0]) */
			for (th1 = MAX_UPARA + 1; th1 >= 0; th1--) {
				th0 = th1;
				min_cost = dpcost[th1] - cbuf_p[th1] + thc_p[0] + thc_p[th0 - tre_p[th0]];
				for (k = 0; k < th1; k++) {
					cost = dpcost[k] - cbuf_p[k] + thc_p[th1 - k] + thc_p[k - tre_p[k]];

					if (cost < min_cost) {
						min_cost = cost;
						th0 = k;
					}
				}

				dpcost[th1] = min_cost + cbuf_p[th1];
				trellis[gr][th1] = th0;

				if (gr == enc->num_group - 1) break;
			}
		}

		th1 = MAX_UPARA + 1;

		for (gr = enc->num_group - 1; gr > 0; gr--) {
			th1 = trellis[gr][th1];
			enc->th[cl][gr - 1] = th1;
		}
	}

	/* set context quantizer */
	for (cl = 0; cl < enc->num_class; cl++) {
		u = 0;

		for (gr = 0; gr < enc->num_group; gr++) {
			for (; u < enc->th[cl][gr]; u++) {
				enc->uquant[cl][u] = gr;
			}
		}
	}

	/* renew groups */
	cost = 0;
	pm_p = enc->pmlist;

	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			cl = enc->class[y][x];
			u = enc->upara[y][x];
			enc->group[y][x] = gr = enc->uquant[cl][u];
			e = enc->encval[y][x];
			prd = enc->prd[y][x];

			if (prd < 0) prd = 0;
			else if (prd > enc->maxprd) prd = enc->maxprd;

			base = enc->bconv[prd];
			pm = pm_p[gr] + enc->fconv[prd];
			cost += pm->cost[base + e] + pm->subcost[base];
		}
	}

	/* optimize probability models */
	if (enc->optimize_loop > 1 && enc->num_pmodel > 1) {
		if (enc->num_pmodel > MAX_UPARA + 2) {
			free(cbuf);
			cbuf = (cost_t **)alloc_2d_array(enc->num_group, enc->num_pmodel, sizeof(cost_t));
		}

		for (gr = 0; gr < enc->num_group; gr++) {
			for (k = 0; k < enc->num_pmodel; k++) {
				cbuf[gr][k] = 0;
			}
		}

		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				gr = enc->group[y][x];
				e = enc->encval[y][x];
				prd = enc->prd[y][x];

				if (prd < 0) prd = 0;
				else if (prd > enc->maxprd) prd = enc->maxprd;

				base = enc->bconv[prd];
				frac = enc->fconv[prd];

				for (k = 0; k < enc->num_pmodel; k++) {
					pm = enc->pmodels[gr][k] + frac;
					cbuf[gr][k] += pm->cost[base + e] + pm->subcost[base];
				}
			}
		}

		for (gr = 0; gr < enc->num_group; gr++) {
			pm = enc->pmodels[gr][0];
			cost = cbuf[gr][0];

			for (k = 1; k < enc->num_pmodel; k++) {
				if (cost > cbuf[gr][k]) {
					cost = cbuf[gr][k];
					pm = enc->pmodels[gr][k];
				}
			}

			pm_p[gr] = pm;
		}

		cost = 0.0;

		for (gr = 0; gr < enc->num_group; gr++) {
			cost += cbuf[gr][pm_p[gr]->id];
		}
	}

	free(cbuf);
	free(dpcost);
	free(trellis);
	return (cost);
}

/*---------------------------- set_prdbuf --------------------------*
 |  Function set_prdbuf
 |
 |  Purpose: Calculates and stores the prediction and residue for all classes for a given block
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN)
 |		prdbuf			--> Struct to store the prediction value (IN/OUT)
 |		errbuf          --> Struct to store the residue value (IN/OUT)
 |		tly				--> Top left corner y position of the current block (IN)
 |		tlx             --> Top left corner x position of the current block (IN)
 |		bufsize			--> Size of the block (IN)
 |
 |  Returns:  ENCODER* --> returns a encoder type structure
 *-------------------------------------------------------------------*/
void set_prdbuf(ENCODER *enc, int **prdbuf, int **errbuf, int tly, int tlx, int bufsize) {
	int x, y, brx, bry, cl, k, prd, *prdbuf_p, *errbuf_p, *coef_p;
	int buf_ptr, org, *org_p, *roff_p;

	//Prediction order determination
	int full_prd_order = enc->prd_order + enc->back_prd_order + enc->for_prd_order;

	//Check the boundaries of the block (mainly for pictures not multiple of bufsize)
	brx = (tlx + bufsize < enc->width) ? (tlx + bufsize) : enc->width;
	bry = (tly + bufsize < enc->height) ? (tly + bufsize) : enc->height;

	//
	for (cl = 0; cl < enc->num_class; cl++) {
		buf_ptr = bufsize * (tly % bufsize) + tlx % bufsize;

		//Run all pixels in the given block
		for (y = tly; y < bry; y++) {
			prdbuf_p = &prdbuf[cl][buf_ptr];
			errbuf_p = &errbuf[cl][buf_ptr];
			buf_ptr += bufsize;
			org_p = &enc->org[1][y][tlx];

			for (x = tlx; x < brx; x++) {
				//If the pixel is of the current class just copy
				if (cl == enc->class[y][x]) {
					*prdbuf_p++ = enc->prd[y][x];
					*errbuf_p++ = enc->err[1][y][x];
					org_p++;
				}
				//If not do the calculations
				else {
					coef_p = enc->predictor[cl];
					roff_p = enc->roff[y][x];
					prd = 0;

					for (k = 0; k < full_prd_order; k++) {
						prd += org_p[*roff_p++] * (*coef_p++);
					}

					org = *org_p++;
					*prdbuf_p++ = prd;

					if (prd < 0) prd = 0;
					else if (prd > enc->maxprd) prd = enc->maxprd;

					prd >>= (enc->coef_precision - 1);
					*errbuf_p++ = enc->econv[org][prd];
				}
			}
		}
	}
}

int find_class(ENCODER *enc, int **prdbuf, int **errbuf, int tly, int tlx, int bry, int brx, int bufsize) {
	cost_t cost, min_cost;
	int x, y, bufptr, cl, min_cl;
	char *class_p;
	int *prd_p, *prdbuf_p, *err_p, *errbuf_p;

	min_cost = 1E8;
	min_cl = 0;

	for (cl = 0; cl < enc->num_class; cl++) {
		bufptr = bufsize * (tly % bufsize) + tlx % bufsize;
		cost = enc->class_cost[enc->mtfbuf[cl]];

		for (y = tly; y < bry; y++) {
			class_p = &enc->class[y][tlx];
			prd_p = &enc->prd[y][tlx];
			prdbuf_p = &prdbuf[cl][bufptr];
			err_p = &enc->err[1][y][tlx];
			errbuf_p = &errbuf[cl][bufptr];
			bufptr += bufsize;

			for (x = tlx; x < brx; x++) {
				*class_p++ = cl;
				*prd_p++ = *prdbuf_p++;
				*err_p++ = *errbuf_p++;
			}
		}

		cost += calc_cost(enc, tly, tlx, bry, brx);

		if (cost < min_cost) {
			min_cost = cost;
			min_cl = cl;
		}
	}

	bufptr = bufsize * (tly % bufsize) + tlx % bufsize;
	for (y = tly; y < bry; y++) {
		class_p = &enc->class[y][tlx];
		prd_p = &enc->prd[y][tlx];
		prdbuf_p = &prdbuf[min_cl][bufptr];
		err_p = &enc->err[1][y][tlx];
		errbuf_p = &errbuf[min_cl][bufptr];
		bufptr += bufsize;
		for (x = tlx; x < brx; x++) {
			*class_p++ = min_cl;
			*prd_p++ = *prdbuf_p++;
			*err_p++ = *errbuf_p++;
		}
	}

	return (min_cl);
}

/*------------------------------ vbs_class --------------------------*
 |  Function vbs_class
 |
 |  Purpose:
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN)
 |		prdbuf			--> Struct to store the prediction value (IN/OUT)
 |		errbuf          --> Struct to store the residue value (IN/OUT)
 |		tly				--> Top left corner y position of the current block (IN)
 |		tlx             --> Top left corner x position of the current block (IN)
 |		blksize			--> Size of the block (IN)
 |		width			--> Image width (IN)
 |		level			--> Level of the quadtree (IN)
 |
 |  Returns:  cost_t --> The class and partition cost
 *-------------------------------------------------------------------*/
cost_t vbs_class(ENCODER *enc, int **prdbuf, int **errbuf, int tly, int tlx, int blksize, int width, int level) {
	int y, x, k, bry, brx, cl, bufsize, bufptr, ctx;
	int mtf_save[MAX_CLASS];
	char **qtmap;
	cost_t cost1, cost2, qtcost;
	char *class_p;
	int *prd_p, *prdbuf_p, *err_p, *errbuf_p;

	//Sets the max block size depending on the current loop
	//and if the quadtree is active
	if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
		bufsize = MAX_BSIZE;
	}
	else {
		bufsize = BASE_BSIZE;
	}

	//Check the boundaries of the block (mainly for pictures not multiple of bufsize)
	brx = (tlx + blksize < enc->width) ? (tlx + blksize) : enc->width;
	bry = (tly + blksize < enc->height) ? (tly + blksize) : enc->height;

	//Safeguard condition
	if (tlx >= brx || tly >= bry) return (0);

	//Auxiliary variable that stores class order
	for (k = 0; k < enc->num_class; k++) {
		mtf_save[k] = enc->mtfbuf[k];
	}

	//In the case of the first loop or if we are not using quadtree
	//Reorders classes to favor the modes chosen by the neighboring blocks.
	mtf_classlabel(enc->class, enc->mtfbuf, tly, tlx, blksize, width, enc->num_class);
	//Finds the class that returns the minimum cost for the block
	cl = find_class(enc, prdbuf, errbuf, tly, tlx, bry, brx, bufsize);
	//Return the class cost
	qtcost = enc->class_cost[enc->mtfbuf[cl]];

	//In the case of the second loop and if we are using quadtree
	//Check the best block partitioning
	if (level > 0) {
		/* context for quad-tree flag */
		ctx = 0;
		qtmap = enc->qtmap[level - 1];
		y = (tly / MIN_BSIZE) >> level;
		x = (tlx / MIN_BSIZE) >> level;

		if (y > 0) {
			if (qtmap[y - 1][x] == 1) ctx++;
			if (brx < width && qtmap[y - 1][x + 1] == 1) ctx++;
		}

		if (x > 0 && qtmap[y][x - 1] == 1) ctx++;

		ctx = ((level - 1) * 4 + ctx) << 1;
		/* Quad-tree partitioning */
		cost1 = calc_cost(enc, tly, tlx, bry, brx) + enc->class_cost[enc->mtfbuf[cl]] + enc->qtflag_cost[ctx];
		blksize >>= 1;

		for (k = 0; k < enc->num_class; k++) {
			enc->mtfbuf[k] = mtf_save[k];
		}

		qtcost = enc->qtflag_cost[ctx + 1];
		qtcost += vbs_class(enc, prdbuf, errbuf, tly, tlx, blksize, width, level - 1);
		qtcost += vbs_class(enc, prdbuf, errbuf, tly, tlx+blksize, blksize, width, level - 1);
		qtcost += vbs_class(enc, prdbuf, errbuf, tly+blksize, tlx, blksize, width, level - 1);
		qtcost += vbs_class(enc, prdbuf, errbuf, tly+blksize, tlx+blksize, blksize, brx, level - 1);
		cost2 = calc_cost(enc, tly, tlx, bry, brx) + qtcost;

		if (cost1 < cost2) {
			blksize <<= 1;

			for (k = 0; k < enc->num_class; k++) {
				enc->mtfbuf[k] = mtf_save[k];
			}

			mtf_classlabel(enc->class, enc->mtfbuf, tly, tlx, blksize, width, enc->num_class);
			qtcost = enc->class_cost[enc->mtfbuf[cl]] + enc->qtflag_cost[ctx];
			bufptr = bufsize * (tly % bufsize) + tlx % bufsize;

			for (y = tly; y < bry; y++) {
				class_p = &enc->class[y][tlx];
				prd_p = &enc->prd[y][tlx];
				prdbuf_p = &prdbuf[cl][bufptr];
				err_p = &enc->err[1][y][tlx];
				errbuf_p = &errbuf[cl][bufptr];
				bufptr += bufsize;

				for (x = tlx; x < brx; x++) {
					*class_p++ = cl;
					*prd_p++ = *prdbuf_p++;
					*err_p++ = *errbuf_p++;
				}
			}

			tly = (tly / MIN_BSIZE) >> level;
			tlx = (tlx / MIN_BSIZE) >> level;
			bry = tly + 1;
			brx = tlx + 1;

			for (; level > 0; level--) {
				qtmap = enc->qtmap[level - 1];
				for (y = tly; y < bry; y++) {
					for (x = tlx; x < brx; x++) {
						qtmap[y][x] = 0;
					}
				}

				tly <<= 1;
				tlx <<= 1;
				bry <<= 1;
				brx <<= 1;
			}
		}
		else {
			qtmap[y][x] = 1;
		}
	}

	return (qtcost);
}

/*--------------------------- optimize_class ------------------------*
 |  Function optimize_class
 |
 |  Purpose: Optimizes the prediction class for each block in the picture
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN/OUT)
 |
 |  Returns:  ENCODER* --> returns a encoder type structure
 *-------------------------------------------------------------------*/
cost_t optimize_class(ENCODER *enc) {
	int y, x, i, blksize, level;
	int **prdbuf, **errbuf;

	//Sets the max block size depending on the current loop
	//and if the quadtree is active
	if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
		level = enc->quadtree_depth;
		blksize = MAX_BSIZE; //32
	}
	else {
		level = 0;
		blksize = BASE_BSIZE; //8
	}

	//Initialize buffer with the class numbers
	for (i = 0; i < enc->num_class; i++) {
		enc->mtfbuf[i] = i;
	}

	//Auxiliary structures, prediction and residue buffer
	prdbuf =(int **) alloc_2d_array(enc->num_class, blksize * blksize, sizeof(int));
	errbuf =(int **) alloc_2d_array(enc->num_class, blksize * blksize, sizeof(int));

	//Cycle to run all picture blocks
	for (y = 0; y < enc->height; y += blksize) {
		for (x = 0; x < enc->width; x += blksize) {
			//Calculates and stores the prediction and residue for all classes for a given block
			set_prdbuf(enc, prdbuf, errbuf, y, x, blksize);
			//
			vbs_class(enc, prdbuf, errbuf, y, x, blksize, enc->width, level);
		}
	}

	//Free auxiliary pointers
	free(errbuf);
	free(prdbuf);

	//Returns cost...
	return (calc_cost(enc, 0, 0, enc->height, enc->width));
}

void optimize_coef(ENCODER *enc, int cl, int pos1, int pos2) {
#define SEARCH_RANGE 11
#define SUBSEARCH_RANGE 3
	cost_t cbuf[SEARCH_RANGE * SUBSEARCH_RANGE], *cbuf_p;
	float *pmcost_p;
	int i, j, k, x, y, df1, df2, base;
	int prd, prd_f, shift, maxprd, *coef_p, *econv_p, *roff_p, *org_p;
	char *class_p;
	img_t *bconv_p, *fconv_p;
	PMODEL *pm, *pm_p;

	cbuf_p = cbuf;
	coef_p = enc->predictor[cl];
	k = 0;

	for (i = 0; i < SEARCH_RANGE; i++) {
		y = coef_p[pos1] + i - (SEARCH_RANGE >> 1);

		if (y < 0) y = -y;

		if (y > MAX_COEF) y = MAX_COEF;

		for (j = 0; j < SUBSEARCH_RANGE; j++) {
			x = coef_p[pos2] - (i - (SEARCH_RANGE >> 1)) - (j - (SUBSEARCH_RANGE >> 1));

			if (x < 0) x = -x;

			if (x > MAX_COEF) x = MAX_COEF;

			cbuf_p[k++] = enc->coef_cost[enc->coef_m[pos1]][y] + enc->coef_cost[enc->coef_m[pos2]][x];
		}
	}

	bconv_p = enc->bconv;
	fconv_p = enc->fconv;
	maxprd = enc->maxprd;
	shift = enc->coef_precision - 1;

	for (y = 0; y < enc->height; y++) {
		class_p = enc->class[y];

		for (x = 0; x < enc->width; x++) {
			if (cl != *class_p++) continue;

			roff_p = enc->roff[y][x];
			prd = enc->prd[y][x];
			org_p = &enc->org[1][y][x];
			pm_p = enc->pmlist[(int)enc->group[y][x]];
			df1 = org_p[roff_p[pos1]];
			df2 = org_p[roff_p[pos2]];
			prd_f = prd - (df1 - df2) * (SEARCH_RANGE >> 1) + df2 * (SUBSEARCH_RANGE >> 1);
			cbuf_p = cbuf;

			if (enc->pm_accuracy < 0) {
				econv_p = enc->econv[*org_p];
				pmcost_p = pm_p->cost;

				for (i = 0; i < SEARCH_RANGE; i++) {
					for (j = 0; j < SUBSEARCH_RANGE; j++) {
						prd = prd_f;

						if (prd < 0) prd = 0;
						else if (prd > maxprd) prd = maxprd;

						(*cbuf_p++) += pmcost_p[econv_p[prd >> shift]];
						prd_f -= df2;
					}

					prd_f += df1 + df2 * (SUBSEARCH_RANGE - 1);
				}
			}
			else {
				for (i = 0; i < SEARCH_RANGE; i++) {
					for (j = 0; j < SUBSEARCH_RANGE; j++) {
						prd = prd_f;

						if (prd < 0) prd = 0;
						else if (prd > maxprd) prd = maxprd;

						base = bconv_p[prd];
						pm = pm_p + fconv_p[prd];
						(*cbuf_p++) += pm->cost[*org_p + base] + pm->subcost[base];
						prd_f -= df2;
					}

					prd_f += df1 + df2 * (SUBSEARCH_RANGE - 1);
				}
			}
		}
	}

	cbuf_p = cbuf;
	j = (SEARCH_RANGE * SUBSEARCH_RANGE) >> 1;

	for (i = 0; i < SEARCH_RANGE * SUBSEARCH_RANGE; i++) {
		if (cbuf_p[i] < cbuf_p[j]) {
			j = i;
		}
	}

	i = (j / SUBSEARCH_RANGE) - (SEARCH_RANGE >> 1);
	j = (j % SUBSEARCH_RANGE) - (SUBSEARCH_RANGE >> 1);
	y = coef_p[pos1] + i;
	x = coef_p[pos2] - i - j;

	if (y < -MAX_COEF) y = -MAX_COEF;
	else if (y > MAX_COEF) y = MAX_COEF;

	if (x < -MAX_COEF) x = -MAX_COEF;
	else if (x > MAX_COEF) x = MAX_COEF;

	i = y - coef_p[pos1];
	j = x - coef_p[pos2];

	if (i != 0 || j != 0) {
		for (y = 0; y < enc->height; y++) {
			class_p = enc->class[y];

			for (x = 0; x < enc->width; x++) {
				if (cl == *class_p++) {
					roff_p = enc->roff[y][x];
					org_p = &enc->org[1][y][x];
					enc->prd[y][x] += org_p[roff_p[pos1]] * i + org_p[roff_p[pos2]] * j;
				}
			}
		}

		coef_p[pos1] += i;
		coef_p[pos2] += j;
	}
}

cost_t optimize_predictor(ENCODER *enc) {
	int cl, k, pos1, pos2;
#ifndef RAND_MAX
#  define RAND_MAX 32767
#endif

	int full_prd_order = enc->prd_order + enc->back_prd_order + enc->for_prd_order;

	for (cl = 0; cl < enc->num_class; cl++) {
		for (k = 0; k < full_prd_order; k++) {
			retry:
			pos1 = (int)(((double)rand() * full_prd_order) / (RAND_MAX+1.0));
			pos2 = (int)(((double)rand() * full_prd_order) / (RAND_MAX+1.0));

			if (pos1 == pos2) goto retry;

			optimize_coef(enc, cl, pos1, pos2);
		}
	}

	predict_region(enc, 0, 0, enc->height, enc->width);
	return (calc_cost(enc, 0, 0, enc->height, enc->width));
}

int putbits(FILE *fp, int n, uint x) {
	static int bitpos = 8;
	static uint bitbuf = 0;
	int bits;

	bits = n;

	if (bits <= 0) return (0);

	while(n >= bitpos) {
		n -= bitpos;

		if (n < 32) {
			bitbuf |= ((x >> n) & (0xff >> (8 - bitpos)));
		}

		putc(bitbuf, fp);
		bitbuf = 0;
		bitpos = 8;
	}

	bitpos -= n;
	bitbuf |= ((x & (0xff >> (8 - n))) << bitpos);

	return (bits);
}

void remove_emptyclass(ENCODER *enc) {
	int cl, i, k, x, y;

	int full_prd_order = enc->prd_order + enc->back_prd_order + enc->for_prd_order;

	for (cl = 0; cl < enc->num_class; cl++) {
		enc->mtfbuf[cl] = 0;
	}

	for (y = 0; y < enc->height; y += MIN_BSIZE) {
		for (x = 0; x < enc->width; x += MIN_BSIZE) {
			cl = enc->class[y][x];
			enc->mtfbuf[cl]++;
		}
	}

	for (i = cl = 0; i < enc->num_class; i++) {
		if (enc->mtfbuf[i] == 0) {
			enc->mtfbuf[i] = -1;
		}
		else {
			enc->mtfbuf[i] = cl++;
		}
	}

	if (cl == enc->num_class) return;	/* no empty class */

	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			i = enc->class[y][x];
			enc->class[y][x] = enc->mtfbuf[i];
		}
	}

	for (i = cl = 0; i < enc->num_class; i++) {
		if (enc->mtfbuf[i] < 0) continue;
		if (cl != i) {
			for (k = 0; k < full_prd_order; k++) {
				enc->predictor[cl][k] = enc->predictor[i][k];
			}
			for (k = 0; k < enc->num_group - 1; k++) {
				enc->th[cl][k] = enc->th[i][k];
			}
		}
		cl++;
	}

	enc->num_class = cl;
}

int write_header(ENCODER *enc, int prd_order[6], int frames, int bframes, int diff, int hevc, FILE *fp) {
	int bits, i;

	bits = putbits(fp, 16, MAGIC_NUMBER);
	bits += putbits(fp, 8, VERSION);
	bits += putbits(fp, 16, enc->width);
	bits += putbits(fp, 16, enc->height);
	bits += putbits(fp, 16, enc->maxval);
	bits += putbits(fp, 16, frames);
	bits += putbits(fp, 8, bframes);
	bits += putbits(fp, 1, hevc);
	bits += putbits(fp, 4, 1);	/* number of components (1 = monochrome) */
	bits += putbits(fp, 6, enc->num_group);
	for (i = 0; i < 6; i++) {
		bits += putbits(fp, 8, prd_order[i]);
	}
	bits += putbits(fp, 1, diff);
	bits += putbits(fp, 8, enc->delta);
	bits += putbits(fp, 6, enc->num_pmodel - 1);
	bits += putbits(fp, 4, enc->coef_precision - 1);
	bits += putbits(fp, 3, enc->pm_accuracy + 1);
	bits += putbits(fp, 1, enc->f_huffman);
	bits += putbits(fp, 1, (enc->quadtree_depth < 0)? 0 : 1);
	bits += putbits(fp, 5, 0);

	return (bits);
}

int write_class(ENCODER *enc, FILE *fp) {
	int bits;

	bits = putbits(fp, 8, enc->num_class);
	return (bits);
}

int encode_golomb(FILE *fp, int m, int v) {
	int bits, p;

	bits = p = (v >> m) + 1;

	while (p > 32) {
		putbits(fp, 32, 0);
		p -= 32;
	}

	putbits(fp, p, 1);	/* prefix code */
	putbits(fp, m, v);

	return (bits + m);
}

void set_qtindex(ENCODER *enc, int *index, uint *hist, int *numidx, int tly, int tlx, int blksize, int width, int level) {
	int i, cl, x, y, ctx;
	char **qtmap;

	if (tly >= enc->height || tlx >= enc->width) return;

	if (level > 0) {
		/* context modeling for quad-tree flag */
		ctx = 0;
		qtmap = enc->qtmap[level - 1];
		y = (tly / MIN_BSIZE) >> level;
		x = (tlx / MIN_BSIZE) >> level;

		if (y > 0) {
			if (qtmap[y - 1][x] == 1) ctx++;
			if (tlx + blksize < width && qtmap[y - 1][x + 1] == 1) ctx++;
		}

		if (x > 0 && qtmap[y][x - 1] == 1) ctx++;

		ctx = ((level - 1) * 4 + ctx) << 1;

		if (qtmap[y][x] == 1) {
			ctx++;
			index[(*numidx)++] = -(ctx + 1);
			enc->qtctx[ctx]++;
			blksize >>= 1;
			set_qtindex(enc, index, hist, numidx, tly, tlx, blksize, width, level - 1);
			set_qtindex(enc, index, hist, numidx, tly, tlx + blksize, blksize, width, level - 1);
			set_qtindex(enc, index, hist, numidx, tly + blksize, tlx, blksize, width, level - 1);
			width = tlx + blksize * 2;

			if (width >= enc->width) width = enc->width;

			set_qtindex(enc, index, hist, numidx, tly + blksize, tlx + blksize, blksize, width, level - 1);

			return;
		}
		else {
			index[(*numidx)++] = -(ctx + 1);
			enc->qtctx[ctx]++;
		}
	}

	cl = enc->class[tly][tlx];
	mtf_classlabel(enc->class, enc->mtfbuf, tly, tlx, blksize, width, enc->num_class);
	i = enc->mtfbuf[cl];
	index[(*numidx)++] = i;
	hist[i]++;

	return;
}

int encode_class(FILE *fp, ENCODER *enc) {
	int i, j, k, numidx, blksize, level, x, y, ctx, bits, *index;
	uint *hist;
	cost_t cost;

#ifndef OPT_SIDEINFO
	if (fp == NULL) return(0);
#endif
	if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
		level = enc->quadtree_depth;
		blksize = MAX_BSIZE;
		numidx = 0;
		y = enc->height / MIN_BSIZE;
		x = enc->width / MIN_BSIZE;

		for (k = 0; k <= level; k++) {
			numidx += y * x;
			y = (y + 1) >> 1;
			x = (x + 1) >> 1;
		}

		for (k = 0; k < QUADTREE_DEPTH << 3; k++) {
			enc->qtctx[k] = 0;
		}

	}
	else {
		level = 0;
		blksize = BASE_BSIZE;
		numidx = enc->height * enc->width / (BASE_BSIZE * BASE_BSIZE);
	}

	hist = (uint *)alloc_mem(enc->num_class * sizeof(uint));
	index = (int *)alloc_mem(numidx * sizeof(int));

	for (i = 0; i < enc->num_class; i++) {
		hist[i] = 0;
		enc->mtfbuf[i] = i;
	}

	numidx = 0;

	for (y = 0; y < enc->height; y += blksize) {
		for (x = 0; x < enc->width; x += blksize) {
			set_qtindex(enc, index, hist, &numidx, y, x, blksize, enc->width, level);
		}
	}

	bits = 0;

	if (enc->f_huffman == 1) {	/* Huffman */
		VLC *vlc;
		vlc = make_vlc(hist, enc->num_class, 16);

		if (fp == NULL) {
			for (i = 0; i < enc->num_class; i++) {
				enc->class_cost[i] = vlc->len[i];
			}

			for (j = 0; j < numidx; j++) {
				i = index[j];

				if (i < 0) {
					bits += 1;
				}
				else {
					bits += enc->class_cost[i];
				}
			}
		}
		else {	/* actually encode */
			for (i = 0; i < enc->num_class; i++) {
				bits += putbits(fp, 4, vlc->len[i] - 1);
			}

			for (j = 0; j < numidx; j++) {
				i = index[j];

				if (i < 0) {
					ctx = -(i + 1);
					putbits(fp, 1, ctx & 1);
				}
				else {
					bits += putbits(fp, vlc->len[i], vlc->code[i]);
				}
			}
		}

		free_vlc(vlc);
	}
	else {			/* Arithmetic */
		PMODEL *pm;
		double p, c;
		int qtree_code[QUADTREE_DEPTH << 2], mtf_code[MAX_CLASS];

		/* context modeling for quad-tree flag */
		if (level > 0) {
			for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
				cost = INT_MAX;

				for (i = k = 0; i < 7; i++) {
					p = qtree_prob[i];
					c = -log(p) * (cost_t)enc->qtctx[(ctx << 1) + 1] - log(1.0 - p) * (cost_t)enc->qtctx[ctx << 1];

					if (c < cost) {
						k = i;
						cost = c;
					}
				}
				p = qtree_prob[qtree_code[ctx] = k];
				enc->qtflag_cost[(ctx << 1) + 1] = -log(p) / log(2.0);
				enc->qtflag_cost[ctx << 1] = -log(1.0 - p) / log(2.0);
			}
		}

		/* quantization of log-transformed probability */
		c = 0.0;

		for (i = 0; i < enc->num_class; i++) {
			c += (double)hist[i];
		}

		for (i = 0; i < enc->num_class; i++) {
			p = (double)hist[i] / c;

			if (p > 0.0) {
				mtf_code[i] = -log(p) / log(2.0) * (PMCLASS_LEVEL / PMCLASS_MAX);

				if (mtf_code[i] >= PMCLASS_LEVEL) {
					mtf_code[i] = PMCLASS_LEVEL - 1;
				}
			}
			else {
				mtf_code[i] = PMCLASS_LEVEL - 1;
			}

			p = exp(-log(2.0) * ((double)mtf_code[i] + 0.5) * PMCLASS_MAX / PMCLASS_LEVEL);
			enc->class_cost[i] = -log(p) / log(2.0);
			hist[i] = p * (1 << 10);

			if (hist[i] <= 0) hist[i] = 1;
		}

		if (fp == NULL) {
			cost = 0.0;

			for (j = 0; j < numidx; j++) {
				i = index[j];

				if (i < 0) {
					ctx = -(i + 1);
					cost += enc->qtflag_cost[ctx];
				}
				else {
					cost += enc->class_cost[i];
				}
			}

			bits = (int)cost;
		}
		else {	/* actually encode */
			PMODEL cpm[1];
			/* additional info. */
			pm = &enc->spm;

			if (level > 0) {
				set_spmodel(pm, 7, -1);

				for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
					i = qtree_code[ctx];
					rc_encode(fp, enc->rc, pm->cumfreq[i], pm->freq[i], pm->cumfreq[pm->size]);
				}
			}

			set_spmodel(pm, PMCLASS_LEVEL, -1);

			for (i = 0; i < enc->num_class; i++) {
				j = mtf_code[i];
				rc_encode(fp, enc->rc, pm->cumfreq[j], pm->freq[j], pm->cumfreq[pm->size]);

				if (pm->cumfreq[pm->size] < MAX_TOTFREQ) {
					for (j = 0; j < pm->size; j++) {
						if (j < mtf_code[i]) {
							pm->freq[j] /= 2;
						}
						else {
							pm->freq[j] *= 2;
						}

						if (pm->freq[j] <= 0) pm->freq[j] = 1;

						pm->cumfreq[j + 1] = pm->cumfreq[j] + pm->freq[j];
					}
				}
			}

			/* set prob. models */
			if (level > 0) {
				pm->size = QUADTREE_DEPTH << 3;

				for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
					i = qtree_code[ctx];
					p = qtree_prob[i];
					pm->freq[(ctx << 1) + 1] = p * (1 << 10);
					p = 1.0 - p;
					pm->freq[(ctx << 1)] = p * (1 << 10);
				}

				for (i = 0; i < pm->size; i++) {
					pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
				}
			}

			cpm->size = enc->num_class;
			cpm->freq = (uint *)alloc_mem((cpm->size * 2 + 1) * sizeof(uint));
			cpm->cumfreq = &(cpm->freq[cpm->size]);
			cpm->cumfreq[0] = 0;

			for (i = 0; i < enc->num_class; i++) {
				cpm->freq[i] = hist[i];
				cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
			}

			for (j = 0; j < numidx; j++) {
				i = index[j];

				if (i < 0) {
					i = -(i + 1);
					ctx = i & (~1);
					rc_encode(fp, enc->rc, pm->cumfreq[i] - pm->cumfreq[ctx], pm->freq[i], pm->cumfreq[ctx + 2] - pm->cumfreq[ctx]);
				}
				else {
					rc_encode(fp, enc->rc, cpm->cumfreq[i], cpm->freq[i], cpm->cumfreq[cpm->size]);
				}
			}

			bits += enc->rc->code;
			enc->rc->code = 0;

			free(cpm->freq);
		}
	}

	free(index);
	free(hist);
	return (bits);
}

int encode_predictor(FILE *fp, ENCODER *enc) {
	int cl, coef, sgn, k, m, min_m, bits = 0;
	cost_t cost, min_cost, t_cost;

	int full_prd_order = enc->prd_order + enc->back_prd_order + enc->for_prd_order;

#ifndef OPT_SIDEINFO
	if (fp == NULL) return(0);
#endif
	t_cost = 0.0;

	for (k = 0; k < full_prd_order; k++) {
		min_cost = INT_MAX;

		for (m = min_m = 0; m < 16; m++) {
			cost = 0.0;

			for (cl = 0; cl < enc->num_class; cl++) {
				coef = enc->predictor[cl][k];

				if (coef < 0) coef = -coef;

				cost += enc->coef_cost[m][coef];
			}

			if (cost < min_cost) {
				min_cost = cost;
				min_m = m;
			}
		}

		t_cost += min_cost;
		enc->coef_m[k] = min_m;
	}

	bits = t_cost;

	if (fp != NULL) {
		bits = 0;

		if (enc->f_huffman == 1) {	/* Huffman */
			for (k = 0; k < full_prd_order; k++) {
				bits += putbits(fp, 4, enc->coef_m[k]);

				for (cl = 0; cl < enc->num_class; cl++) {
					coef = enc->predictor[cl][k];
					sgn = (coef < 0)? 1 : 0;

					if (coef < 0) coef = -coef;

					bits += encode_golomb(fp, enc->coef_m[k], coef);

					if (coef != 0) {
						bits += putbits(fp, 1, sgn);
					}
				}
			}
		} else {			/* Arithmetic */
			PMODEL *pm;
			pm = &enc->spm;

			for (k = 0; k < full_prd_order; k++) {
				set_spmodel(pm, MAX_COEF + 1, enc->coef_m[k]);
				rc_encode(fp, enc->rc, enc->coef_m[k], 1, 16);

				for (cl = 0; cl < enc->num_class; cl++) {
					coef = enc->predictor[cl][k];
					sgn = (coef < 0)? 1 : 0;

					if (coef < 0) coef = -coef;

					rc_encode(fp, enc->rc, pm->cumfreq[coef],  pm->freq[coef], pm->cumfreq[pm->size]);

					if (coef > 0) {
						rc_encode(fp, enc->rc, sgn, 1, 2);
					}
				}
			}

			bits = enc->rc->code;
			enc->rc->code = 0;
		}
	}

	return (bits);
}

int encode_threshold(FILE *fp, ENCODER *enc) {
	int cl, gr, i, k, m, min_m, bits = 0;
	cost_t cost, min_cost;
	PMODEL *pm;

#ifndef OPT_SIDEINFO
	if (fp == NULL) return(0);
#endif
	if (enc->f_huffman == 1) {	/* Huffman */
		min_cost = INT_MAX;

		for (m = min_m = 0; m < 16; m++) {
			bits = 0;

			for (cl = 0; cl < enc->num_class; cl++) {
				k = 0;

				for (gr = 1; gr < enc->num_group; gr++) {
					i = enc->th[cl][gr - 1] - k;
					bits++;

					if (i > 0) {
						bits += ((i - 1) >> m) + m + 1;
					}

					k += i;

					if (k > MAX_UPARA) break;
				}
			}

			if ((cost = bits) < min_cost) {
				min_cost = cost;
				min_m = m;
			}
		}

		if (fp == NULL) {
			enc->th_cost[0] = 1.0;

			for (i = 1; i < MAX_UPARA + 2; i++) {
				enc->th_cost[i] =((i - 1) >> min_m) + min_m + 1 + 1;
			}

			bits = min_cost;
		}
		else {
			bits = putbits(fp, 4, min_m);

			for (cl = 0; cl < enc->num_class; cl++) {
				k = 0;

				for (gr = 1; gr < enc->num_group; gr++) {
					i = enc->th[cl][gr - 1] - k;

					if (i == 0) {
						bits += putbits(fp, 1, 0);
					} else {
						bits += putbits(fp, 1, 1);
						bits += encode_golomb(fp, min_m, i - 1);
					}

					k += i;

					if (k > MAX_UPARA) break;
				}
			}

			if (enc->num_pmodel > 1) {
				for (k = 1; (1 << k) < enc->num_pmodel; k++);

				for (gr = 0; gr < enc->num_group; gr++) {
					pm = enc->pmlist[gr];
					bits += putbits(fp, k, pm->id);
				}
			}
		}
	}
	else {			/* Arithmetic */
		double p;
		pm = &enc->spm;
		min_cost = INT_MAX;

		for (m = min_m = 0; m < 16; m++) {
			set_spmodel(pm, MAX_UPARA + 2, m);
			cost = 0.0;

			for (cl = 0; cl < enc->num_class; cl++) {
				k = 0;

				for (gr = 1; gr < enc->num_group; gr++) {
					i = enc->th[cl][gr - 1] - k;
					p = (double)pm->freq[i] / (pm->cumfreq[pm->size - k]);
					cost += -log(p);
					k += i;

					if (k > MAX_UPARA) break;
				}
			}

			cost /= log(2.0);

			if (cost < min_cost) {
				min_cost = cost;
				min_m = m;
			}
		}

		set_spmodel(pm, MAX_UPARA + 2, min_m);
		p = log(pm->cumfreq[MAX_UPARA + 2]);

		if (fp == NULL) {
			for (i = 0; i < MAX_UPARA + 2; i++) {
				enc->th_cost[i] = (p - log(pm->freq[i])) / log(2.0);
			}

			bits = min_cost;
		}
		else {
			rc_encode(fp, enc->rc, min_m, 1, 16);

			for (cl = 0; cl < enc->num_class; cl++) {
				k = 0;

				for (gr = 1; gr < enc->num_group; gr++) {
					i = enc->th[cl][gr - 1] - k;
					rc_encode(fp, enc->rc, pm->cumfreq[i],  pm->freq[i], pm->cumfreq[pm->size - k]);
					k += i;

					if (k > MAX_UPARA) break;
				}
			}

			if (enc->num_pmodel > 1) {
				for (gr = 0; gr < enc->num_group; gr++) {
					pm = enc->pmlist[gr];
					rc_encode(fp, enc->rc, pm->id, 1, enc->num_pmodel);
				}
			}

			bits = enc->rc->code;
			enc->rc->code = 0;
		}
	}

	return (bits);
}

void print_contexts(int modo, uint cumfreq, uint freq, uint totfreq, int max, int min, int prd) {
	if (modo == 0) system("rm /tmp/encoder_context.txt");

	FILE *fp;
	fp = fileopen("/tmp/encoder_context.txt", "a");

	if (modo != -1) {
		fprintf(fp, "Frame: %d\n", modo);
		fprintf(fp, "Cum. Freq.\tFreq.\tTot. Freq.\tMax\tMin\tPrd:\n");
	}
	else {
		fprintf(fp, "%d\t\t%d\t%d\t\t%d\t%d\t%d\n", cumfreq, freq, totfreq, max, min, prd);
		//fprintf(fp, "%d\n", prd);
	}

	fclose(fp);
}

int encode_image(FILE *fp, ENCODER *enc) {
	int x, y, e, prd, base, bits, gr, cumbase;
	PMODEL *pm;

	bits = 0;
	if (enc->f_huffman == 1) {	/* Huffman */
		VLC *vlc;

		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				gr = enc->group[y][x];
				e = enc->encval[y][x];
				pm = enc->pmlist[gr];
				vlc = &enc->vlcs[gr][pm->id];
				bits += putbits(fp, vlc->len[e], vlc->code[e]);
			}
		}
	}
	else {			/* Arithmetic */
		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				gr = enc->group[y][x];
				prd = enc->prd[y][x];

				if (prd < 0) prd = 0;
				else if (prd > enc->maxprd) prd = enc->maxprd;

				e = enc->encval[y][x];
				base = enc->bconv[prd];
				pm = enc->pmlist[gr] + enc->fconv[prd];
				cumbase = pm->cumfreq[base];

//				print_contexts(-1, 0, 0, pm->cumfreq[base + enc->maxval + 1] - cumbase, base + enc->maxval + 1, base, prd);

				rc_encode(fp, enc->rc, pm->cumfreq[base + e] - cumbase, pm->freq[base + e], pm->cumfreq[base + enc->maxval + 1] - cumbase);
			}
		}

		rc_finishenc(fp, enc->rc);
		bits += enc->rc->code;
	}

	return (bits);
}

/* 	remove_ext: removes the "extension" from a file spec.
   	   mystr is the string to process.
   	   dot is the extension separator.
   	   sep is the path separator (0 means to ignore).
	Returns an allocated string identical to the original but
   	   with the extension removed. It must be freed when you're
   	   finished with it.
	If you pass in NULL or the new string can't be allocated,
   	   it returns NULL. */
char *remove_ext (char* mystr, char dot, char sep) {
	char *retstr, *lastdot, *lastsep;

	// Error checks and allocate string.
	if (mystr == NULL)
		return NULL;
	if ((retstr = malloc (strlen (mystr) + 1)) == NULL)
		return NULL;

	// Make a copy and find the relevant characters.
	strcpy (retstr, mystr);
	lastdot = strrchr (retstr, dot);
	lastsep = (sep == 0) ? NULL : strrchr (retstr, sep);

	// If it has an extension separator.
	if (lastdot != NULL) {
		// and it's before the extension separator.
		if (lastsep != NULL) {
			if (lastsep < lastdot) {
				// then remove it.
				*lastdot = '\0';
			}
		}
		else {
			// Has extension separator with no path separator.
			*lastdot = '\0';
		}
	}

	// Return the modified string.
	return retstr;
}

void print_results(FILE *res, int frames, int height, int width, int header, int *class_info, int *predictors, int *thresholds, int *errors, int *extrainfo, int prate, int pframes) {
	int f, rate = 0, total_bits = 0;

	printf("\n---------------------------------\n");
	printf("Header info.\t :%10d bits\n", header);
	fprintf(res, "Header info.\t :%10d bits\n", header);
	fprintf(res, "Frame\t\tBits\t\tBpp\n");

	for (f = 0; f < frames; f++) {
		// Results output
		printf("---------------------------------\n");
		printf("frame [%03d]\n", f);
		printf("class info\t :%10d bits\n", class_info[f]);
		printf("predictors\t :%10d bits\n", predictors[f]);
		printf("thresholds\t :%10d bits\n", thresholds[f]);
		printf("pred. errors\t :%10d bits\n", errors[f]);

		rate = class_info[f] + predictors[f] + thresholds[f] + errors[f];

		if (extrainfo != NULL) {
			printf("extra info\t :%10d bits\n", extrainfo[f]);
			rate += extrainfo[f];
		}

		total_bits += rate;

		printf("---------------------------------\n");
		printf("total frame [%2d] :%10d bits\n", f, rate);
		printf("frame coding rate:%10.5f b/p\n", (double)rate / (height * width));
		fprintf(res, "%d\t%10d\t%10.5f\n", f, rate, (double)rate / (height * width));
	}

	total_bits += header;

	printf("---------------------------------\n");
	printf("total\t\t :%10d bits\n", total_bits);
	printf("total coding rate:%10.5f b/p\n", (double)total_bits / (height * width * frames));

	if (pframes != 0) {
		printf("\nAverage I bitrate:%10.5f bits\n", (double) (class_info[0] + predictors[0] + thresholds[0] + errors[0]) / (height * width));
		printf("Average P bitrate:%10.5f bits (%d frames)\n", (double) prate / (pframes * height * width), pframes);
		printf("Average B bitrate:%10.5f bits (%d frames)\n", (double) (total_bits - (class_info[0] + predictors[0] + thresholds[0] + errors[0]) - prate) / ((frames - pframes - 1) * height * width), frames - pframes - 1);
	}

	fprintf(res, "------------------------------------------------\n");
	fprintf(res, "Total:\n\t%10d\t%10.5f\n", total_bits, (double)total_bits / (height * width * frames));

	if (pframes != 0) {
		fprintf(res, "\nAverage I bitrate:%10.5f bits\n", (double) (class_info[0] + predictors[0] + thresholds[0] + errors[0]) / (height * width));
		fprintf(res, "Average P bitrate:%10.5f bits (%d frames)\n", (double) prate / (pframes * height * width), pframes);
		fprintf(res, "Average B bitrate:%10.5f bits (%d frames)\n", (double) (total_bits - (class_info[0] + predictors[0] + thresholds[0] + errors[0]) - prate) / ((frames - pframes - 1) * height * width), frames - pframes - 1);
	}

	if (extrainfo != NULL) {
		total_bits = 0;
		for (f = 0; f < frames; f++) {
			total_bits += extrainfo[f];
		}

		printf("\nExtra info due to pixelwise difference: %d bits --> %f bpp.", total_bits, (double)total_bits / (height * width * frames));
	}
}

IMAGE* calc_diff(IMAGE* ref, IMAGE* cur, char **extra_info, int *num_pels) {
	int x, y;
	// Image allocation
	IMAGE *diff = alloc_image(ref->width, ref->height, 255);
	int aux, conta = 0;
	char *aux_extra_info;

	aux_extra_info = (char *) alloc_mem(ref->width * ref->height * sizeof(char));

	(*extra_info) = NULL;

	for (y = 0; y < ref->height; y++) {
		for (x = 0; x < ref->width; x++) {
			aux = cur->val[y][x] - ref->val[y][x] + 127;

			if (aux > 0 && aux < 255) {
				diff->val[y][x] = aux;
			}
			else if (aux <= 0) {
				diff->val[y][x] = 0;

				aux_extra_info[conta] = aux;

				conta++;
			}
			else if (aux >= 255) {
				diff->val[y][x] = 255;

				aux_extra_info[conta] = aux - 255;

				conta++;
			}
		}
	}

	(*extra_info) = (char *) alloc_mem(conta * sizeof(char));

	for (y = 0; y < conta; y++) {
		(*extra_info)[y] = aux_extra_info[y];
	}

	*num_pels = conta;

	free(aux_extra_info);
	free(ref->val);
	free(ref);
	free(cur->val);
	free(cur);

	return diff;
}

int encode_extra_info(FILE *fp, char *extra_info, int num_pels) {
	int bits = 0, i;

	bits = putbits(fp, 16, num_pels);

	for (i = 0; i < num_pels; i++) {
		bits += putbits(fp, 8, extra_info[i]);
	}

	return bits;
}

void print_predictors(ENCODER *enc, int f) {
	int y, x;

	if (f == 0) system("rm encoder_predictors.txt");

	FILE *teste;
	teste = fileopen("encoder_predictors.txt", "a");

	//fprintf(teste, "Frame: %d\n\n", f);
	for (y = 0; y < enc->num_class; y++) {
		//fprintf(teste, "intra: ");
		for (x = 0; x < enc->prd_order; x++) {
			fprintf(teste, "%d ", enc->predictor[y][x]);
		}
		//fprintf(teste, "|| inter: ");
		for (x = 0; x < enc->back_prd_order; x++) {
			fprintf(teste, "%d ", enc->predictor[y][enc->prd_order + x]);
		}fprintf(teste, "\n");
	}fprintf(teste, "\n");

	fclose(teste);
}

int main(int argc, char **argv) {
	// Variable declaration
	cost_t cost, min_cost, side_cost;
	int i, j, f, k, x, y, cl, **prd_save, **th_save, **error = NULL;
	char **class_save;
	IMAGE *video[3] = {NULL, NULL, NULL};
	ENCODER *enc;
	double elapse = 0.0;
	int f_mmse = 0;
	int f_optpred = 0;
	int f_huffman = 0;
	int quadtree_depth = QUADTREE_DEPTH;
	int num_class = NUM_CLASS;
	int num_group = NUM_GROUP;
	int bframes = 0;
	int prd_order[6] = {INTRA_PRD_ORDER, INTRA_PRD_ORDER, INTER_PRD_ORDER, INTRA_PRD_ORDER, INTER_PRD_ORDER, INTER_PRD_ORDER};
	int coef_precision = COEF_PRECISION;
	int num_pmodel = NUM_PMODEL;
	int pm_accuracy = PM_ACCURACY;
	int max_iteration = MAX_ITERATION;
	int frames = FRAMES;
	int height = 0, width = 0;
	char *infile, *outfile;
	char resfile[100];
	FILE *fp, *res;
	//Print results variables
	int header = 0, *class_info, *predictors, *thresholds, *errors, *extrainfo = NULL;
	int delta = 1;
	int hevc = 1;
	int diff = 0, num_pels = 0;
	char *extra_info = NULL;

	cpu_time();
	setbuf(stdout, 0);
	infile = outfile = NULL;

	// Read input parameters
	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
			case 'H':
				height = atoi(argv[++i]);
				if (height <= 0) {
					exit(-1);
				}
				break;
			case 'W':
				width = atoi(argv[++i]);
				if (width <= 0) {
					exit(-2);
				}
				break;
			case 'F':
				frames = atoi(argv[++i]);
				if (frames <= 0) {
					frames = FRAMES;
				}
				break;
			case 'p':
				diff = 1;
				break;
			case 'M':
				num_class = atoi(argv[++i]);
				if (num_class <= 0 || num_class > 63) {
					num_class = NUM_CLASS;
				}
				break;
			case 'G':
				bframes = atoi(argv[++i]);
				if (bframes <= 0) {
					bframes = 0;
				}
				else {
					bframes = bframes + 1;
				}
				break;
			case 'K':
				for (j = 0; j < 6; j++) {
					prd_order[j] = atoi(argv[++i]);

					if (prd_order[j] < 0 || prd_order[j] > 72) {
						if (j == 0 || j == 1 || j == 3) {
							prd_order[j] = INTRA_PRD_ORDER;
						}
						else {
							prd_order[j] = INTER_PRD_ORDER;
						}
					}
				}
				break;
			case 'P':
				coef_precision = atoi(argv[++i]);
				if (coef_precision <= 0 || coef_precision > 16) {
					coef_precision = COEF_PRECISION;
				}
				break;
			case 'V':
				num_pmodel = atoi(argv[++i]);
				if (num_pmodel <= 0 || num_pmodel > 64) {
					num_pmodel = NUM_PMODEL;
				}
				break;
			case 'A':
				pm_accuracy = atoi(argv[++i]);
				if (pm_accuracy < -1 || pm_accuracy > 6) {
					pm_accuracy = PM_ACCURACY;
				}
				break;
			case 'I':
				max_iteration = atoi(argv[++i]);
				if (max_iteration <= 0) {
					max_iteration = MAX_ITERATION;
				}
				break;
			case 'm':
				f_mmse = 1;
				break;
			case 'o':
				f_optpred = 1;
				break;
			case 'h':
				f_huffman = 1;
				break;
			case 'f':
				quadtree_depth = -1;
				break;
			case 'D':
				delta = atoi(argv[++i]);
				if (delta <= 0) {
					delta = 1;
				}
				break;
			case 'B':
				hevc = atoi(argv[++i]);
				if (hevc != 0 && hevc != 1) {
					hevc = 0;
				}
				break;
			default:
				fprintf(stderr, "Unknown option: %s!\n", argv[i]);
				exit (1);
			}
		}
		else {
			if (infile == NULL) {
				infile = argv[i];
			}
			else {
				outfile = argv[i];
			}
		}
	}

	//If the Huffman coding is used "turns off" the accuracy of the probabilities model
	if (f_huffman == 1) pm_accuracy = -1;
	//The accuracy of the probabilities model can't be greater than the coefficients precision
	if (pm_accuracy > coef_precision) pm_accuracy = coef_precision;

	// Help and version
	if (infile == NULL || outfile == NULL) {
		printf(BANNER"\n", 0.1 * VERSION);
		printf("usage: encmrp [options] infile outfile\n");
		printf("options:\n");
		printf("    -H num  Height*\n");
		printf("    -W num  Width*\n");
		printf("    -F num  Frames*\n");
		printf("    -p 		Pixel difference prediction\n");
		printf("    -M num  Number of predictors [%d]\n", num_class);
		printf("    -G num  Number of frames between references (number of B frames) [%d]\n", bframes);
		printf("    -K num  Prediction order [%d %d %d %d %d %d]\n", prd_order[0], prd_order[1], prd_order[2], prd_order[3], prd_order[4], prd_order[5]);
		printf("    -P num  Precision of prediction coefficients (fractional bits) [%d]\n", coef_precision);
		printf("    -V num  Number of probability models [%d]\n", num_pmodel);
		printf("    -A num  Accuracy of probability models [%d]\n", pm_accuracy);
		printf("    -I num  Maximum number of iterations [%d]\n", max_iteration);
		printf("    -m      Use MMSE predictors\n");
		printf("    -h      Use Huffman coding\n");
		printf("    -f      Fixed block-size for adaptive prediction\n");
		printf("    -o      Further optimization of predictors (experimental)\n");
		printf("infile:     Input file (must be in a raw PGM format)\n");
		printf("outfile:    Output file\n");
		printf("\nNote: * stands for a mandatory option.\n");
		exit(0);
	}

	if (bframes != 0 && frames < bframes + 1) {
		printf("Error: not enough frames to use biprediction, the program will now exit.\n");
		exit(-2);
	}
	
	// Open output file
	fp = fileopen(outfile, "wb");

	// If the number of classes was not defined it is now
	k = width * height;
	if (num_class < 0) {
		num_class = 10.4E-5 * k + 13.8;

		if (num_class > MAX_CLASS) {
			num_class = MAX_CLASS;
		}
	}

	// If the prediction order was not defined it is now
	for (j = 0; j < 6; j++) {
		if (prd_order[j] < 0) {
			prd_order[j] = 12.0E-5 * k + 17.2;

			for (i = 1; i < 8; i++) {
				if (prd_order[j] < (i + 1) * (i + 1)) {
					prd_order[j] = i * (i + 1);
					break;
				}
			}

			if (i >= 8) prd_order[j] = 72;
		}
	}

	printf("\nMRP-Video Encoder\n\n");
	// Print file characteristics to screen
	printf("%s -> %s (%dx%dx%d)\n", infile, outfile, width, height, frames);
	// Print coding parameters to screen
	printf("M = %d, P = %d, V = %d, A = %d, D = %d, p = %s\n\n", num_class, coef_precision, num_pmodel, pm_accuracy, delta, (diff == 1) ? "on": "off");
	printf("Number of B frames: %d\nPrediction order:\n\tFrame I: %d\n\tFrame P: %d %d\n\tFrame B: %d %d %d\n\n", bframes == 0 ? 0 : bframes - 1, prd_order[0], prd_order[1], prd_order[2], prd_order[3], prd_order[4], prd_order[5]);

	//Allocation of print results variables
	errors 	   = (int *) alloc_mem(frames * sizeof(int));
	class_info = (int *) alloc_mem(frames * sizeof(int));
	predictors = (int *) alloc_mem(frames * sizeof(int));
	thresholds = (int *) alloc_mem(frames * sizeof(int));
	extrainfo  = (int *) alloc_mem(frames * sizeof(int));
	for (f = 0; f < frames; f++) {
		errors[f] = 0;
		class_info[f] = 0;
		predictors[f] = 0;
		thresholds[f] = 0;
		extrainfo[f] = 0;
	}

	char *aux = remove_ext(outfile, '.', '/');
	sprintf(resfile, "res_%s.txt", aux);
	free(aux);

	res = fileopen(resfile, "w");
	fprintf(res, "MRP-Video version %d encoding results\n", VERSION);
	fprintf(res, "\tEncoded file: %s\n", infile);
	fprintf(res, "\tDimensions: %dx%dx%d\n", width, height, frames);
	fprintf(res, "---------------------------------------------\n");

	//Loop to encode all the frames
	if (bframes == 0) {
		for (f = 0; f < frames; f++) {
			// Creates new ENCODER structure
			if (f == 0) {
				// Read input file
				video[1] = read_yuv(infile, height, width, f);

				enc = init_encoder(video[1], NULL, NULL, NULL, NULL, num_class, num_group, prd_order[0], 0, 0, coef_precision, f_huffman, quadtree_depth, num_pmodel, pm_accuracy, delta);
			}
			else {
				// Read input file
				if (diff == 0) {
					video[0] = video[1];
					video[1] = read_yuv(infile, height, width, f);
				}
				else {
					video[0] = video[1];
					video[1] = calc_diff(read_yuv(infile, height, width, f - 1), read_yuv(infile, height, width, f), &extra_info, &num_pels);
				}

				enc = init_encoder(video[1], video[0], NULL, error, NULL, num_class, num_group, prd_order[1], prd_order[2], 0, coef_precision, f_huffman, quadtree_depth, num_pmodel, pm_accuracy, delta);

				free(error);
			}

			int full_prd_order = enc->prd_order + enc->back_prd_order;

			//Initiation of the reference offset
			enc->roff = init_ref_offset(video[1], enc->prd_order, enc->back_prd_order, 0);

			// Initialize probability models
			enc->pmodels = init_pmodels(enc->num_group, enc->num_pmodel, enc->pm_accuracy, NULL, enc->sigma, enc->maxval + 1);

			// Huffman coding
			if (enc->f_huffman == 1) {
				enc->vlcs = init_vlcs(enc->pmodels, enc->num_group, enc->num_pmodel);
			}

			// Set cost model
			set_cost_model(enc, f_mmse);

			// Set cost model
			init_class(enc);

			//Auxiliary variables
			prd_save = (int **)alloc_2d_array(enc->num_class, (full_prd_order), sizeof(int));
			th_save = (int **)alloc_2d_array(enc->num_class, enc->num_group, sizeof(int));
			class_save = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));

			/* 1st loop */
			//Loop type
			enc->optimize_loop = 1;
			min_cost = INT_MAX;

			for (i = j = 0; i < max_iteration; i++) {
				printf("[%2d] cost =", i);
				cost = design_predictor(enc, f_mmse);
				printf(" %d ->", (int)cost);
				cost = optimize_group(enc);
				printf(" %d ->", (int)cost);
				cost = optimize_class(enc);
				printf(" %d", (int)cost);

				if (cost < min_cost) {
					 printf(" *\n");
					min_cost = cost;
					j = i;

					for (y = 0; y < enc->height; y++) {
						for (x = 0; x < enc->width; x++) {
							class_save[y][x] = enc->class[y][x];
						}
					}

					for (cl = 0; cl < enc->num_class; cl++) {
						for (k = 0; k < full_prd_order; k++) {
							prd_save[cl][k] = enc->predictor[cl][k];
						}

						for (k = 0; k < enc->num_group; k++) {
							th_save[cl][k] = enc->th[cl][k];
						}
					}
				}
				else {
					printf("\n");
				}

				if (i - j >= EXTRA_ITERATION) break;
				elapse += cpu_time();
			}

			//Restore values
			for (y = 0; y < enc->height; y++) {
				for (x = 0; x < enc->width; x++) {
					enc->class[y][x] = class_save[y][x];
				}
			}
			for (cl = 0; cl < enc->num_class; cl++) {
				for (k = 0; k < full_prd_order; k++) {
					enc->predictor[cl][k] = prd_save[cl][k];
				}
				for (k = 0; k < enc->num_group; k++) {
					enc->th[cl][k] = th_save[cl][k];
				}
			}

			set_cost_rate(enc);
			predict_region(enc, 0, 0, enc->height, enc->width);
			cost = calc_cost(enc, 0, 0, enc->height, enc->width);

			printf("Frame: %d\n\t1st optimization --> Cost: %d\n", f, (int)cost);

			/* 2nd loop */
			//Loop type
			enc->optimize_loop = 2;
			min_cost = INT_MAX;

			for (i = j = 0; i < max_iteration; i++) {
				if (f_optpred) {
					cost = optimize_predictor(enc);
				}

				side_cost = encode_predictor(NULL, enc);
				cost = optimize_group(enc);
				side_cost += encode_threshold(NULL, enc);
				cost = optimize_class(enc);
				side_cost += encode_class(NULL, enc);
				cost += side_cost;

				if (cost < min_cost) {
					min_cost = cost;
					j = i;

					if (f_optpred) {
						for (y = 0; y < enc->height; y++) {
							for (x = 0; x < enc->width; x++) {
								class_save[y][x] = enc->class[y][x];
							}
						}
						for (cl = 0; cl < enc->num_class; cl++) {
							for (k = 0; k < full_prd_order; k++) {
								prd_save[cl][k] = enc->predictor[cl][k];
							}
							for (k = 0; k < enc->num_group; k++) {
								th_save[cl][k] = enc->th[cl][k];
							}
						}
					}
				}

				if (f_optpred) {
					if (i - j >= EXTRA_ITERATION) break;
				}
				else {
					if (i > j) break;
				}
				elapse += cpu_time();
			}

			if (f_optpred) {
				for (y = 0; y < enc->height; y++) {
					for (x = 0; x < enc->width; x++) {
						enc->class[y][x] = class_save[y][x];
					}
				}

				for (cl = 0; cl < enc->num_class; cl++) {
					for (k = 0; k < full_prd_order; k++) {
						enc->predictor[cl][k] = prd_save[cl][k];
					}
					i = 0;

					for (k = 0; k < enc->num_group; k++) {
						enc->th[cl][k] = th_save[cl][k];
						for (; i < enc->th[cl][k]; i++) {
							enc->uquant[cl][i] = k;
						}
					}
				}

				predict_region(enc, 0, 0, enc->height, enc->width);
				calc_cost(enc, 0, 0, enc->height, enc->width);
				optimize_class(enc);
			}

			remove_emptyclass(enc);

			printf("\t2nd optimization --> Cost: %d (%d)", (int)cost, (int)side_cost);
			printf(" --> M = %d\n", enc->num_class);

			if (f == 0) {
				header = write_header(enc, prd_order, frames, bframes, diff, 0, fp);
			}

			class_info[f] = write_class(enc, fp);

			if (enc->f_huffman == 0) {
				enc->rc = rc_init();
			}

			class_info[f] += encode_class(fp, enc);
			predictors[f] = encode_predictor(fp, enc);
			thresholds[f] = encode_threshold(fp, enc);
			errors[f] = encode_image(fp, enc);

			if (f > 0 && diff == 1) {
				extrainfo[f] += encode_extra_info(fp, extra_info, num_pels);
			}

			if (f > 0) {
				free(video[0]->val);
				free(video[0]);
			}

			free(class_save);
			free(prd_save);
			free(th_save);

			error = get_enc_err(enc, 1);

			free_encoder(enc);
		}

		free(video[1]->val);
		free(video[1]);

		if (f_huffman == 1) {
			putbits(fp, 7, 0);	/* flush remaining bits */
		}

		fclose(fp);

		if (diff == 1) {
			print_results(res, frames, height, width, header, class_info, predictors, thresholds, errors, extrainfo, 0, 0);
		}
		else {
			print_results(res, frames, height, width, header, class_info, predictors, thresholds, errors, NULL, 0, 0);
		}

		free(error);
	}
	else {
		int **back_ref_error, **for_ref_error;
		int back_reference, for_reference, first_frame = 0, conta = 0, final;
		int ***keep_error;
		int (*bref)[5];
		int prate = 0, pframes = 0;

		if (hevc == 0 || bframes == 2) {
			hevc = 0;
			back_reference = 0;
			for_reference = bframes;
		}
		else if (hevc == 1) {
			keep_error = (int ***)alloc_mem((bframes + 1) * sizeof(int **));

			for (f = 0; f < bframes + 1; f++) keep_error[f] = NULL;

			bref = select_bref(bframes);

			final = ((frames - 1) % bframes != 0) ? (bframes * ((int)((frames - 1) / bframes)) + 1) : frames;
		}

		f = 0;

		while(f < frames) {
			// Creates new ENCODER structure
			if (f == 0) {
				// Read input file
				video[1] = read_yuv(infile, height, width, f);

				enc = init_encoder(video[1], NULL, NULL, NULL, NULL, num_class, num_group, prd_order[0], 0, 0, coef_precision, f_huffman, quadtree_depth, num_pmodel, pm_accuracy, delta);
			}
			else if ((f - first_frame) % bframes == 0) {
				enc = init_encoder(video[1], video[0], NULL, back_ref_error, NULL, num_class, num_group, prd_order[1], prd_order[2], 0, coef_precision, f_huffman, quadtree_depth, num_pmodel, pm_accuracy, delta);
			}
			else {
				enc = init_encoder(video[1], video[0], video[2], back_ref_error, for_ref_error, num_class, num_group, prd_order[3], prd_order[4], prd_order[5], coef_precision, f_huffman, quadtree_depth, num_pmodel, pm_accuracy, delta);
			}

			int full_prd_order = enc->prd_order + enc->back_prd_order + enc->for_prd_order;

			//Initiation of the reference offset
			enc->roff = init_ref_offset(video[1], enc->prd_order, enc->back_prd_order, enc->for_prd_order);

			// Initialize probability models
			enc->pmodels = init_pmodels(enc->num_group, enc->num_pmodel, enc->pm_accuracy, NULL, enc->sigma, enc->maxval + 1);

			// Huffman coding
			if (enc->f_huffman == 1) {
				enc->vlcs = init_vlcs(enc->pmodels, enc->num_group, enc->num_pmodel);
			}

			// Set cost model
			set_cost_model(enc, f_mmse);

			// Set cost model
			init_class(enc);

			//Auxiliary variables
			prd_save = (int **)alloc_2d_array(enc->num_class, (full_prd_order), sizeof(int));
			th_save = (int **)alloc_2d_array(enc->num_class, enc->num_group, sizeof(int));
			class_save = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));

			/* 1st loop */
			//Loop type
			enc->optimize_loop = 1;
			min_cost = INT_MAX;

			for (i = j = 0; i < max_iteration; i++) {
				cost = design_predictor(enc, f_mmse);
				cost = optimize_group(enc);
				cost = optimize_class(enc);

				if (cost < min_cost) {
					min_cost = cost;
					j = i;

					for (y = 0; y < enc->height; y++) {
						for (x = 0; x < enc->width; x++) {
							class_save[y][x] = enc->class[y][x];
						}
					}

					for (cl = 0; cl < enc->num_class; cl++) {
						for (k = 0; k < full_prd_order; k++) {
							prd_save[cl][k] = enc->predictor[cl][k];
						}

						for (k = 0; k < enc->num_group; k++) {
							th_save[cl][k] = enc->th[cl][k];
						}
					}

				}

				if (i - j >= EXTRA_ITERATION) break;
				elapse += cpu_time();
			}

			//Restore values
			for (y = 0; y < enc->height; y++) {
				for (x = 0; x < enc->width; x++) {
					enc->class[y][x] = class_save[y][x];
				}
			}
			for (cl = 0; cl < enc->num_class; cl++) {
				for (k = 0; k < full_prd_order; k++) {
					enc->predictor[cl][k] = prd_save[cl][k];
				}
				for (k = 0; k < enc->num_group; k++) {
					enc->th[cl][k] = th_save[cl][k];
				}
			}

			set_cost_rate(enc);
			predict_region(enc, 0, 0, enc->height, enc->width);
			cost = calc_cost(enc, 0, 0, enc->height, enc->width);

			printf("Frame: %d\n\t1st optimization --> Cost: %d\n", f, (int)cost);

			/* 2nd loop */
			//Loop type
			enc->optimize_loop = 2;
			min_cost = INT_MAX;

			for (i = j = 0; i < max_iteration; i++) {
				if (f_optpred) {
					cost = optimize_predictor(enc);
				}

				side_cost = encode_predictor(NULL, enc);
				cost = optimize_group(enc);
				side_cost += encode_threshold(NULL, enc);
				cost = optimize_class(enc);
				side_cost += encode_class(NULL, enc);
				cost += side_cost;

				if (cost < min_cost) {
					min_cost = cost;
					j = i;

					if (f_optpred) {
						for (y = 0; y < enc->height; y++) {
							for (x = 0; x < enc->width; x++) {
								class_save[y][x] = enc->class[y][x];
							}
						}
						for (cl = 0; cl < enc->num_class; cl++) {
							for (k = 0; k < full_prd_order; k++) {
								prd_save[cl][k] = enc->predictor[cl][k];
							}
							for (k = 0; k < enc->num_group; k++) {
								th_save[cl][k] = enc->th[cl][k];
							}
						}
					}
				}

				if (f_optpred) {
					if (i - j >= EXTRA_ITERATION) break;
				}
				else {
					if (i > j) break;
				}
				elapse += cpu_time();
			}

			if (f_optpred) {
				for (y = 0; y < enc->height; y++) {
					for (x = 0; x < enc->width; x++) {
						enc->class[y][x] = class_save[y][x];
					}
				}

				for (cl = 0; cl < enc->num_class; cl++) {
					for (k = 0; k < full_prd_order; k++) {
						enc->predictor[cl][k] = prd_save[cl][k];
					}
					i = 0;

					for (k = 0; k < enc->num_group; k++) {
						enc->th[cl][k] = th_save[cl][k];
						for (; i < enc->th[cl][k]; i++) {
							enc->uquant[cl][i] = k;
						}
					}
				}

				predict_region(enc, 0, 0, enc->height, enc->width);
				calc_cost(enc, 0, 0, enc->height, enc->width);
				optimize_class(enc);
			}

			remove_emptyclass(enc);

			printf("\t2nd optimization --> Cost: %d (%d)", (int)cost, (int)side_cost);
			printf(" --> M = %d\n", enc->num_class);

			if (f == 0) {
				header = write_header(enc, prd_order, frames, bframes, diff, hevc, fp);
			}

			class_info[f] = write_class(enc, fp);

			if (enc->f_huffman == 0) {
				enc->rc = rc_init();
			}

			class_info[f] += encode_class(fp, enc);
			predictors[f] = encode_predictor(fp, enc);
			thresholds[f] = encode_threshold(fp, enc);
			errors[f] = encode_image(fp, enc);

			if (f > 0 && diff == 1) {
				extrainfo[f] += encode_extra_info(fp, extra_info, num_pels);
			}

			if ((((f - first_frame) % bframes == 0) | (f == frames - 1)) && f != 0) {
				prate += class_info[f] + predictors[f] + thresholds[f] + errors[f];
				pframes++;
			}

			free(class_save);
			free(prd_save);
			free(th_save);

			//Next frame selection
			if (hevc == 0) {
				if (f == 0) {
					f = bframes;

					back_ref_error = get_enc_err(enc, 1);

					video[0] = copy_yuv(video[1]);
					safefree_yuv(&video[1]);

					if (diff == 0) {
						video[1] = read_yuv(infile, height, width, f);
					}
					else {
						video[1] = calc_diff(read_yuv(infile, height, width, f - 1), read_yuv(infile, height, width, f), &extra_info, &num_pels);
					}
				}
				else if (f == for_reference) {
					if (f == back_reference + 1) {
						f = frames;
					}
					else {
						f = back_reference + 1;

						for_ref_error = get_enc_err(enc, 1);

						if (video[2] != NULL) safefree_yuv(&video[2]);

						video[2] = copy_yuv(video[1]);
						safefree_yuv(&video[1]);

						if (diff == 0) {
							video[1] = read_yuv(infile, height, width, f);
						}
						else {
							video[1] = calc_diff(read_yuv(infile, height, width, f - 1), read_yuv(infile, height, width, f), &extra_info, &num_pels);
						}
					}
				}
				else if (back_reference < f && f < for_reference - 1) {
					f++;

					safefree_yuv(&video[1]);

					if (diff == 0) {
						video[1] = read_yuv(infile, height, width, f);
					}
					else {
						video[1] = calc_diff(read_yuv(infile, height, width, f - 1), read_yuv(infile, height, width, f), &extra_info, &num_pels);
					}
				}
				else if (f == for_reference - 1) {
					if (f + 1 == frames - 1) {
						f = frames;
					}
					else if (f + bframes + 1 < frames) {
						f = f + bframes + 1;

						safefree((void **)&back_ref_error);
						back_ref_error = for_ref_error;
						back_reference = for_reference;
						for_reference = f;
					}
					else {
						f = frames - 1;

						safefree((void **)&back_ref_error);
						back_ref_error = for_ref_error;
						back_reference = for_reference;
						for_reference = f;
					}

					safefree_yuv(&video[0]);
					video[0] = copy_yuv(video[2]);

					safefree_yuv(&video[1]);

					if (diff == 0) {
						video[1] = read_yuv(infile, height, width, f);
					}
					else {
						video[1] = calc_diff(read_yuv(infile, height, width, f - 1), read_yuv(infile, height, width, f), &extra_info, &num_pels);
					}
				}
			}
			else if (hevc == 1) {
				if (f == 0) {
					f = bref[conta][0];

					keep_error[0] = get_enc_err(enc, 1);
					for_ref_error = NULL;
					video[2] = NULL;
					back_ref_error = keep_error[0];

					safefree_yuv(&video[1]);
					if (diff == 0) {
						video[0] = read_yuv(infile, height, width, f + bref[conta][1]);
						video[1] = read_yuv(infile, height, width, f);
					}
					else {
						video[0] = read_yuv(infile, height, width, f + bref[conta][1]);
						video[1] = calc_diff(read_yuv(infile, height, width, f - 1), read_yuv(infile, height, width, f), &extra_info, &num_pels);
					}

					conta++;
				}
				else {
					if (conta + 1 <= bframes) {
						keep_error[f - first_frame] = get_enc_err(enc, 1);

						if (f + bref[conta][0] > 2 * first_frame && first_frame != 0) {
							f = first_frame + bref[conta][0];
						}
						else {
							f = f + bref[conta][0];
						}

						safefree_yuv(&video[1]);
						safefree_yuv(&video[0]);
						if (diff == 0) {
							video[0] = read_yuv(infile, height, width, f + bref[conta][1]);
							video[1] = read_yuv(infile, height, width, f);
						}
						else {
							if (f + bref[conta][1] == 0) {
								video[0] = read_yuv(infile, height, width, f + bref[conta][1]);
							}
							else {
								video[0] = calc_diff(read_yuv(infile, height, width, f + bref[conta][1] - 1), read_yuv(infile, height, width, f + bref[conta][1]), &extra_info, &num_pels);
							}
							video[1] = calc_diff(read_yuv(infile, height, width, f - 1), read_yuv(infile, height, width, f), &extra_info, &num_pels);
						}

						if (conta < bframes + 1) back_ref_error = keep_error[bref[conta][3]];

						if (bref[conta][2] != -1) {
							if (video[2] != NULL) safefree_yuv(&video[2]);
							if (diff == 0) {
								video[2] = read_yuv(infile, height, width, f + bref[conta][2]);
							}
							else {
								video[2] = calc_diff(read_yuv(infile, height, width, f + bref[conta][2] - 1), read_yuv(infile, height, width, f + bref[conta][2]), &extra_info, &num_pels);
							}

							if (conta < bframes + 1) for_ref_error = keep_error[bref[conta][4]];
						}
						else {
							if (video[2] != NULL) safefree_yuv(&video[2]);
							for_ref_error = NULL;
						}

						conta++;
					}
					else if (conta + 1 > bframes) {
						if (first_frame + 2 * bframes >= frames) {
							if (f + 1 < final && final != frames) {
								int aux_bframes = bframes;
								bframes = frames - final;

								conta = 1;

								safefree_yuv(&video[0]);
								safefree_yuv(&video[1]);
								if (diff == 0) {
									video[0] = read_yuv(infile, height, width, f + 1);
									video[1] = read_yuv(infile, height, width, frames - 1);
								}
								else {
									video[0] = calc_diff(read_yuv(infile, height, width, f), read_yuv(infile, height, width, f + 1), &extra_info, &num_pels);
									video[1] = calc_diff(read_yuv(infile, height, width, frames - 2), read_yuv(infile, height, width, frames - 1), &extra_info, &num_pels);
								}

								first_frame = f + 1;
								f = frames - 1;
								final = frames;

								if (bframes == 1) {
									conta = aux_bframes;

									safefree_yuv(&video[0]);
									safefree_yuv(&video[1]);
									if (diff == 0) {
										video[0] = read_yuv(infile, height, width, frames - 2);
										video[1] = read_yuv(infile, height, width, frames - 1);
									}
									else {
										video[0] = calc_diff(read_yuv(infile, height, width, frames - 3), read_yuv(infile, height, width, frames - 2), &extra_info, &num_pels);
										video[1] = calc_diff(read_yuv(infile, height, width, frames - 2), read_yuv(infile, height, width, frames - 1), &extra_info, &num_pels);
									}

									first_frame = frames - 2;
								}
								else {
									bref = select_bref(bframes);
								}

								safefree_yuv(&video[2]);
								video[2] = NULL;
								for_ref_error = NULL;

								back_ref_error = keep_error[aux_bframes];
								safefree((void **)&keep_error[0]);
								keep_error[0] = keep_error[aux_bframes];

								for (i = 1; i < aux_bframes; i++) {
									safefree((void **)&keep_error[i]);
								}
							}
							else {
								break;
							}
						}
						else {
							conta = 1;
							f = first_frame + 2 * bframes;

							safefree_yuv(&video[0]);
							safefree_yuv(&video[1]);
							if (diff == 0) {
								video[0] = read_yuv(infile, height, width, f + bref[0][1]);
								video[1] = read_yuv(infile, height, width, f);
							}
							else {
								video[0] = calc_diff(read_yuv(infile, height, width, f + bref[0][1] - 1), read_yuv(infile, height, width, f + bref[0][1]), &extra_info, &num_pels);
								video[1] = calc_diff(read_yuv(infile, height, width, f - 1), read_yuv(infile, height, width, f), &extra_info, &num_pels);
							}

							safefree_yuv(&video[2]);
							video[2] = NULL;
							for_ref_error = NULL;

							back_ref_error = keep_error[bframes];
							safefree((void **)&keep_error[0]);
							keep_error[0] = keep_error[bframes];

							for (i = 1; i < bframes; i++) {
								safefree((void **)&keep_error[i]);
							}

							first_frame = first_frame + bframes;
						}
					}
				}
			}

			free_encoder(enc);
			enc = NULL;
		}

		if (enc != NULL) free_encoder(enc);

		if (video[1] != video[0] && video[1] != video[2]) {
			if (video[1] != NULL) safefree_yuv(&video[1]);
		}
		if (video[0] != video[2]) {
			safefree_yuv(&video[0]);
		}
		if (video[2] != NULL) safefree_yuv(&video[2]);

		if (f_huffman == 1) {
			putbits(fp, 7, 0);	/* flush remaining bits */
		}

		if (hevc == 0 || bframes == 2) {
			if (back_ref_error != for_ref_error) {
				safefree((void **)&back_ref_error);
			}
			safefree((void **)&for_ref_error);
		}
		else if (hevc == 1) {
			for (i = 0; i < bframes + 1; i++) {
				safefree((void **)&keep_error[i]);
			}

			safefree((void **)&keep_error);
		}

		fclose(fp);

		if (diff == 1) {
			print_results(res, frames, height, width, header, class_info, predictors, thresholds, errors, extrainfo, prate, pframes);
		}
		else {
			print_results(res, frames, height, width, header, class_info, predictors, thresholds, errors, NULL, prate, pframes);
		}
	}

	free(class_info);
	free(predictors);
	free(thresholds);
	free(errors);
	free(extrainfo);

	elapse += cpu_time();

	fprintf(res, "\nCPU time: %.2f sec.\n\n", elapse);
	fclose(res);
	printf("\ncpu time: %.2f sec.\n", elapse);

	return (0);
}
