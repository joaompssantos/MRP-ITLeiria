#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "mrp.h"
#include "mrp_config.h"

extern POINT dyx[];
extern POINT idyx[];
extern POINT tridyx[];
extern double sigma_h[], sigma_a[];
extern double qtree_prob[];

uint getbits(FILE *fp, int n) {
	static int bitpos = 0;
	static uint bitbuf = 0;
	int x = 0;

	if (n <= 0) return (0);

	while (n > bitpos) {
		n -= bitpos;
		x = (x << bitpos) | bitbuf;
		bitbuf = (uint) (getc(fp) & 0xff);
		bitpos = 8;
	}

	bitpos -= n;
	x = (x << n) | (bitbuf >> bitpos);
	bitbuf &= ((1 << bitpos) - 1);

	return (uint) x;
}

int read_class(FILE *fp) {
	return(getbits(fp, 8));
}

DECODER *init_decoder(FILE *fp, int **back_ref_error, int **for_ref_error, int width, int height, int maxval, int num_comp, int num_group, int prd_order, int back_prd_order, int for_prd_order, int mi_size, int *mi_prd_order, int num_pmodel, int coef_precision, int pm_accuracy, int f_huffman, int quadtree_depth, int delta, int depth) {
	DECODER *dec;
	int i;

	dec = (DECODER *)alloc_mem(sizeof(DECODER));

	dec->width = width;
	dec->height = height;
	dec->maxval = maxval;
	dec->num_comp = num_comp;
	dec->num_group = num_group;

	dec->prd_order = prd_order;
	dec->back_prd_order = back_prd_order;
	dec->for_prd_order = for_prd_order;

	dec->mi_prd_order[UP] = mi_prd_order[UP];
	dec->mi_prd_order[LEFT] = mi_prd_order[LEFT];
	dec->mi_prd_order[RDIAG] = mi_prd_order[RDIAG];
	dec->mi_prd_order[LDIAG] = mi_prd_order[LDIAG];

	dec->full_prd_order = dec->prd_order + dec->back_prd_order + dec->for_prd_order + dec->mi_prd_order[UP] +
						  dec->mi_prd_order[LEFT] + dec->mi_prd_order[RDIAG] + dec->mi_prd_order[LDIAG];

	dec->mi_size = mi_size;
	dec->num_pmodel = num_pmodel;
	dec->coef_precision = coef_precision;
	dec->pm_accuracy = pm_accuracy;
	dec->f_huffman = f_huffman;
	dec->quadtree_depth = quadtree_depth;
	dec->maxprd = dec->maxval << dec->coef_precision;
	dec->delta = delta;
	dec->depth = depth;

	dec->num_class = read_class(fp);

	dec->predictor = (int **)alloc_2d_array(dec->num_class, dec->full_prd_order, sizeof(int));

	dec->err = (int ***)alloc_3d_array(dec->height, dec->width, 3, sizeof(int));

	if (back_ref_error != NULL) {
		int y, x;

		for (y = 0; y < dec->height; y++) {
			for (x = 0; x < dec->width; x++) {
				dec->err[0][y][x] = back_ref_error[y][x];
			}
		}
	}

	if (for_ref_error != NULL) {
		int y, x;

		for (y = 0; y < dec->height; y++) {
			for (x = 0; x < dec->width; x++) {
				dec->err[2][y][x] = for_ref_error[y][x];
			}
		}
	}

	//Initiation of the reference offset
	dec->intra_roff = init_ref_offset(dec->height, dec->width, INTRA_PRED, dec->prd_order, 0);
	// Keeps the weights used for the residue encoding context
	dec->intra_ctx_weight = init_ctx_weight(INTRA_PRED, dec->prd_order, dec->delta);

	if (dec->back_prd_order > 0) {
		dec->back_roff = init_ref_offset(dec->height, dec->width, BACK_PRED, dec->back_prd_order, 0);
		dec->back_ctx_weight = init_ctx_weight(BACK_PRED, dec->back_prd_order, dec->delta);
	}

	if (dec->for_prd_order > 0) {
		dec->for_roff = init_ref_offset(dec->height, dec->width, FOR_PRED, dec->for_prd_order, 0);
		dec->for_ctx_weight = init_ctx_weight(FOR_PRED, dec->prd_order, dec->delta);
	}

	if (dec->mi_prd_order[UP] > 0) {
		dec->mi_up_roff = init_ref_offset(dec->height, dec->width, MI_UP_PRED, dec->mi_prd_order[UP], dec->mi_size);
		dec->mi_up_ctx_weight = init_ctx_weight(MI_UP_PRED, dec->mi_prd_order[UP], dec->delta);

		dec->mi_borders_ctx_weight = init_ctx_weight(MI_BORDERS_PRED, dec->mi_prd_order[UP], dec->delta);
	}

	if (dec->mi_prd_order[LEFT] > 0) {
		dec->mi_left_roff = init_ref_offset(dec->height, dec->width, MI_LEFT_PRED, dec->mi_prd_order[LEFT], dec->mi_size);
		dec->mi_left_ctx_weight = init_ctx_weight(MI_LEFT_PRED, dec->mi_prd_order[LEFT], dec->delta);

		if (dec->mi_prd_order[LEFT] > dec->mi_prd_order[UP]) {
			dec->mi_borders_ctx_weight = init_ctx_weight(MI_BORDERS_PRED, dec->mi_prd_order[LEFT], dec->delta);
		}
	}

	if (dec->mi_prd_order[LDIAG] > 0) {
		dec->mi_ldiag_roff = init_ref_offset(dec->height, dec->width, MI_LDIAG_PRED, dec->mi_prd_order[LDIAG], dec->mi_size);
		dec->mi_ldiag_ctx_weight = init_ctx_weight(MI_LDIAG_PRED, dec->mi_prd_order[LDIAG], dec->delta);

		if (dec->mi_prd_order[LDIAG] > dec->mi_prd_order[UP] && dec->mi_prd_order[LDIAG] > dec->mi_prd_order[LEFT]) {
			dec->mi_borders_ctx_weight = init_ctx_weight(MI_BORDERS_PRED, dec->mi_prd_order[LDIAG], dec->delta);
		}
	}

	if (dec->mi_prd_order[RDIAG] > 0) {
		dec->mi_rdiag_roff = init_ref_offset(dec->height, dec->width, MI_RDIAG_PRED, dec->mi_prd_order[RDIAG], dec->mi_size);
		dec->mi_rdiag_ctx_weight = init_ctx_weight(MI_RDIAG_PRED, dec->mi_prd_order[RDIAG], dec->delta);

		if (dec->mi_prd_order[RDIAG] > dec->mi_prd_order[UP] && dec->mi_prd_order[RDIAG] > dec->mi_prd_order[LEFT] &&
			dec->mi_prd_order[RDIAG] > dec->mi_prd_order[LDIAG]) {
			dec->mi_borders_ctx_weight = init_ctx_weight(MI_BORDERS_PRED, dec->mi_prd_order[RDIAG], dec->delta);
		}
	}

	if (dec->f_huffman == 0) {
		dec->rc = rc_init();
		rc_startdec(fp, dec->rc);
	}

	// Quadtree map
	if (dec->quadtree_depth > 0) {
		int x, y, xx, yy;

		yy = (dec->height + MAX_BSIZE - 1) / MAX_BSIZE;
		xx = (dec->width + MAX_BSIZE - 1) / MAX_BSIZE;

		for (i = dec->quadtree_depth - 1; i >= 0; i--) {
			dec->qtmap[i] = (char **)alloc_2d_array(yy, xx, sizeof(char));

			for (y = 0; y < yy; y++) {
				for (x = 0; x < xx; x++) {
					dec->qtmap[i][y][x] = 0;
				}
			}

			yy <<= 1;
			xx <<= 1;
		}
	}

	// Class and uquant arrays
	dec->class = (char **)alloc_2d_array(dec->height, dec->width, sizeof(char));

	dec->uquant = (char **)alloc_2d_array(dec->num_class, MAX_UPARA + 1, sizeof(char));

	if (dec->num_pmodel > 1) {
		dec->pm_idx = (int *)alloc_mem(dec->num_group * sizeof(int));
	}
	else {
		dec->pm_idx = NULL;
	}

	dec->spm.freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
	dec->spm.cumfreq = &(dec->spm.freq[MAX_SYMBOL]);

	//Huffman coding
	if (dec->f_huffman == 1) {
		dec->sigma = sigma_h;
	}
	else {
		dec->sigma = sigma_a;
	}

	dec->mtfbuf = (int *)alloc_mem(dec->num_class * sizeof(int));

	return (dec);
}

int **get_dec_err(DECODER *dec, int pos) {
	int y, x;
	int **error = (int **)alloc_2d_array(dec->height + 1, dec->width, sizeof(int));

	for (y = 0; y < dec->height; y++) {
		for (x = 0; x < dec->width; x++) {
			error[y][x] = dec->err[pos][y][x];
		}
	}

	return error;
}


void free_decoder(DECODER *dec) {
	int i, j, gr, num_subpm;

	free(dec->predictor);
	free(dec->err);

	free(dec->intra_ctx_weight);
    free_ref_offset(dec->height, dec->width, INTRA_PRED, dec->prd_order, 0, dec->intra_roff);

	if(dec->back_prd_order > 0) {
	    free(dec->back_ctx_weight);
        free_ref_offset(dec->height, dec->width, BACK_PRED, dec->back_prd_order, 0, dec->back_roff);
	}

	if(dec->for_prd_order > 0) {
	    free(dec->for_ctx_weight);
        free_ref_offset(dec->height, dec->width, FOR_PRED, dec->for_prd_order, 0, dec->for_roff);

    }

	if(dec->mi_prd_order[UP] > 0) {
	    free(dec->mi_up_ctx_weight);
        free_ref_offset(dec->height, dec->width, MI_UP_PRED, dec->mi_prd_order[UP], dec->mi_size, dec->mi_up_roff);
	}

	if(dec->mi_prd_order[LEFT] > 0) {
	    free(dec->mi_left_ctx_weight);
        free_ref_offset(dec->height, dec->width, MI_LEFT_PRED, dec->mi_prd_order[LEFT], dec->mi_size, dec->mi_left_roff);
	}

	if(dec->mi_prd_order[LDIAG] > 0) {
	    free(dec->mi_ldiag_ctx_weight);
        free_ref_offset(dec->height, dec->width, MI_LDIAG_PRED, dec->mi_prd_order[LDIAG], dec->mi_size, dec->mi_ldiag_roff);
	}

	if(dec->mi_prd_order[RDIAG] > 0) {
	    free(dec->mi_rdiag_ctx_weight);
        free_ref_offset(dec->height, dec->width, MI_RDIAG_PRED, dec->mi_prd_order[RDIAG], dec->mi_size, dec->mi_rdiag_roff);
	}

	if(dec->mi_prd_order[UP] > 0 || dec->mi_prd_order[LEFT] > 0 || dec->mi_prd_order[LDIAG] > 0 || dec->mi_prd_order[RDIAG] > 0) free(dec->mi_borders_ctx_weight);


    if (dec->quadtree_depth > 0) {
		for (i = dec->quadtree_depth - 1; i >= 0; i--) {
			free(dec->qtmap[i]);
		}
	}

	free(dec->class);
	free(dec->uquant);

	if (dec->num_pmodel > 1) {
		free(dec->pm_idx);
	}

	free(dec->spm.freq);

	if (dec->pm_accuracy < 0) {
		num_subpm = 1;
	}
	else {
		num_subpm = 1 << dec->pm_accuracy;
	}
	for (gr = 0; gr < dec->num_group; gr++) {
		for (j = 0; j < num_subpm; j++) {
			PMODEL *aux = &dec->pmodels[gr][0][j];
			free(aux->freq);
		}
	}
	free(dec->pmodels[0][0]);
	free(dec->pmodels);

	free(dec->mtfbuf);

	if (dec->f_huffman == 0) {
		free(dec->rc);
	}

	free(dec);
}

void read_header(FILE *fp, int *width, int *height, int *maxval, int *frames, int *depth, int *bframes, int *num_comp,
				 int *num_group, int prd_order[6], int *mi_size, int mi_prd_order[4], int *num_pmodel,
				 int *coef_precision, int *pm_accuracy, int *f_huffman, int *quadtree_depth, int *delta, int *hevc,
				 int *hist_bytes) {
	int i = 0;

	if (getbits(fp, 16) != MAGIC_NUMBER) {
		fprintf(stderr, "Not a compressed file!\n");
		exit(1);
	}

	*width = getbits(fp, 16);
	*height = getbits(fp, 16);
	*maxval = getbits(fp, 16);
	*frames = getbits(fp, 16);
	*depth = getbits(fp, 5);
	*bframes = getbits(fp, 8);
	*hevc = getbits(fp, 1);
	*num_comp = getbits(fp, 4);
	*num_group = getbits(fp, 6);

	for (i = 0; i < 6; i++) {
		prd_order[i] = getbits(fp, 8);
	}

	*mi_size = getbits(fp, 8);
	for (i = 0; i < 4; i++) {
		mi_prd_order[i] = getbits(fp, 8);
	}

	*delta = getbits(fp, 8);
	*num_pmodel = getbits(fp, 6) + 1;
	*coef_precision = getbits(fp, 4) + 1;
	*pm_accuracy = getbits(fp, 4) - 1;
	*f_huffman = getbits(fp, 1);
	*quadtree_depth = (getbits(fp, 1))? QUADTREE_DEPTH : -1;
	*hist_bytes = getbits(fp, 16);
}

int decode_vlc(FILE *fp, VLC *vlc) {
	int i, k, min, off;
	uint code;

	code = 0;
    min = off = 0;

	for (i = 0; i < vlc->max_len; i++) {
		code = (code << 1) | getbits(fp, 1);
		k = vlc->off[i];

		if (k < 0) {
			min <<= 1;
		}
		else {
			if (code <= vlc->code[vlc->index[k]]) break;

			min = (vlc->code[vlc->index[k]] + 1) << 1;
			off = k + 1;
		}
	}

	i = off + code - min;

	return (vlc->index[i]);
}

int decode_golomb(FILE *fp, int m) {
	int v = 0;

	while (getbits(fp, 1) == 0) {
		v++;
	}

	v = (v << m) | getbits(fp, m);

	return (v);
}

void decode_predictor(FILE *fp, DECODER *dec) {
	int k, m, cl, coef, sgn;

	if (dec->f_huffman == 1) {
		for (k = 0; k < dec->full_prd_order; k++) {
			m = getbits(fp, 4);

			for (cl = 0; cl < dec->num_class; cl++) {
				coef = decode_golomb(fp, m);

				if (coef > 0) {
					sgn = getbits(fp, 1);

					if (sgn) {
						coef = -coef;
					}
				}

				dec->predictor[cl][k] = coef;
			}
		}
	}
	else {
		PMODEL *pm;

		pm = &dec->spm;
		pm->size = MAX_COEF + 18;
		pm->cumfreq[MAX_COEF + 2] = 0;

		for (k = MAX_COEF + 2; k < pm->size; k++) {
			pm->freq[k] = 1;
			pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
		}

		for (k = 0; k < dec->full_prd_order; k++) {
			m = rc_decode(fp, dec->rc, pm, MAX_COEF + 2, MAX_COEF + 18) - (MAX_COEF + 2);
			set_spmodel(pm, MAX_COEF + 1, m);

			for (cl = 0; cl < dec->num_class; cl++) {
				coef = rc_decode(fp, dec->rc, pm, 0, MAX_COEF + 1);

				if (coef > 0) {
					sgn = rc_decode(fp, dec->rc, pm, MAX_COEF+2, MAX_COEF+4) - (MAX_COEF + 2);

					if (sgn) {
						coef = -coef;
					}
				}

				dec->predictor[cl][k] = coef;
			}
		}
	}
}

void decode_threshold(FILE *fp, DECODER *dec) {
	int cl, gr, m, u, k;

	if (dec->f_huffman == 1) {
		m = getbits(fp, 4);

		for (cl = 0; cl < dec->num_class; cl++) {
			k = u = 0;

			for (gr = 0; gr < dec->num_group; gr++) {
				if (k > MAX_UPARA || gr == dec->num_group - 1) {
					k = MAX_UPARA + 1;
				}
				else {
					if (getbits(fp, 1)) k += decode_golomb(fp, m) + 1;
				}

				for (; u < k; u++) dec->uquant[cl][u] = (char) gr;
			}
		}

		if (dec->num_pmodel > 1) {
			for (k = 1; (1 << k) < dec->num_pmodel; k++);

			for (gr = 0; gr < dec->num_group; gr++) {
				dec->pm_idx[gr] = getbits(fp, k);
			}
		}
	}
	else {
		PMODEL *pm;

		pm = &dec->spm;
		pm->size = 16;
		pm->cumfreq[0] = 0;

		for (k = 0; k < pm->size; k++) {
			pm->freq[k] = 1;
			pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
		}

		m = rc_decode(fp, dec->rc, pm, 0, pm->size);
		set_spmodel(pm, MAX_UPARA + 2, m);

		for (cl = 0; cl < dec->num_class; cl++) {
			k = u = 0;

			for (gr = 0; gr < dec->num_group; gr++) {
				if (k > MAX_UPARA || gr == dec->num_group - 1) {
					k = MAX_UPARA + 1;
				}
				else {
					k += rc_decode(fp, dec->rc, pm, 0, pm->size - k);
				}

				for (; u < k; u++) dec->uquant[cl][u] = (char) gr;
			}
		}

		if (dec->num_pmodel > 1) {
			pm->size = dec->num_pmodel;
			pm->freq[0] = 0;

			for (k = 0; k < pm->size; k++) {
				pm->freq[k] = 1;
				pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
			}

			for (gr = 0; gr < dec->num_group; gr++) {
				dec->pm_idx[gr] = rc_decode(fp, dec->rc, pm, 0, pm->size);
			}
		}
	}
}

void decode_qtindex(FILE *fp, DECODER *dec, VLC *vlc, PMODEL *cpm, int tly, int tlx, int blksize, int width, int level) {
	int i, cl, y, x, bry, brx, ctx;
	char **qtmap;
	PMODEL *pm;

	brx = (tlx + blksize < dec->width) ? (tlx + blksize) : dec->width;
	bry = (tly + blksize < dec->height) ? (tly + blksize) : dec->height;

	if (tlx >= brx || tly >= bry) return;

	if (level > 0) {
		ctx = 0;
		qtmap = dec->qtmap[level - 1];
		y = (tly / MIN_BSIZE) >> level;
		x = (tlx / MIN_BSIZE) >> level;

		if (y > 0) {
			if (qtmap[y - 1][x] == 1) ctx++;
			if (brx < width && qtmap[y - 1][x + 1] == 1) ctx++;
		}

		if (x > 0 && qtmap[y][x - 1] == 1) ctx++;

		ctx = ((level - 1) * 4 + ctx) << 1;

		if (dec->f_huffman == 1) {
			i = getbits(fp, 1);
		}
		else {
			pm = &dec->spm;
			i = rc_decode(fp, dec->rc, pm, ctx, ctx + 2) - ctx;
		}

		if (i == 1) {
			qtmap[y][x] = 1;
			blksize >>= 1;
			decode_qtindex(fp, dec, vlc, cpm, tly, tlx, blksize, width, level - 1);
			decode_qtindex(fp, dec, vlc, cpm, tly, tlx + blksize, blksize, width, level - 1);
			decode_qtindex(fp, dec, vlc, cpm, tly + blksize, tlx, blksize, width, level - 1);
			decode_qtindex(fp, dec, vlc, cpm, tly + blksize, tlx + blksize, blksize, brx, level - 1);

			return;
		}
	}

	if (dec->f_huffman == 1) {
		i = decode_vlc(fp, vlc);
	}
	else {
		i = rc_decode(fp, dec->rc, cpm, 0, cpm->size);
	}

	mtf_classlabel(dec->class, dec->mtfbuf, tly, tlx,blksize, width, dec->num_class);

	for (cl = 0; cl < dec->num_class; cl++) {
		if (dec->mtfbuf[cl] == i) break;
	}

	for (y = tly; y < bry; y++) {
		for (x = tlx; x < brx; x++) {
			dec->class[y][x] = (char) cl;
		}
	}
}

void decode_class(FILE *fp, DECODER *dec) {
	int i, j, x, y, blksize, level;
	VLC *vlc;
	PMODEL *pm, cpm[1];

	if (dec->quadtree_depth >= 0) {
		level = dec->quadtree_depth;
		blksize = MAX_BSIZE;
	}
	else {
		level = 0;
		blksize = BASE_BSIZE;
	}

	if (dec->f_huffman == 1) {
		vlc = (VLC *)alloc_mem(sizeof(VLC));
		vlc->size = dec->num_class;
		vlc->max_len = 16;
		vlc->len = (int *)alloc_mem(vlc->size * sizeof(int));

		for (i = 0; i < vlc->size; i++) {
			vlc->len[i] = getbits(fp, 4) + 1;
		}

		gen_huffcode(vlc);
	}
	else {
		double p;
		int ctx;
		int qtree_code[QUADTREE_DEPTH << 2], mtf_code[MAX_CLASS];

		pm = &dec->spm;
		if (dec->quadtree_depth > 0) {
			set_spmodel(pm, 7, -1);

			for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
				qtree_code[ctx] = rc_decode(fp, dec->rc, pm, 0, pm->size);
			}
		}

		set_spmodel(pm, PMCLASS_LEVEL, -1);

		for (i = 0; i < dec->num_class; i++) {
			mtf_code[i] = rc_decode(fp, dec->rc, pm, 0, pm->size);

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
				pm->freq[(ctx << 1) + 1] = (uint) (p * (1 << 10));
				p = 1.0 - p;
				pm->freq[(ctx << 1)] = (uint) (p * (1 << 10));
			}

			for (i = 0; i < pm->size; i++) {
				pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
			}
		}

		cpm->size = dec->num_class;
		cpm->freq = (uint *)alloc_mem((cpm->size * 2 + 1) * sizeof(uint));
		cpm->cumfreq = &cpm->freq[cpm->size];
		cpm->cumfreq[0] = 0;

		for (i = 0; i < dec->num_class; i++) {
			p = exp(-log(2.0) * ((double)mtf_code[i] + 0.5) * PMCLASS_MAX/PMCLASS_LEVEL);
			cpm->freq[i] = (uint) (p * (1 << 10));

			if (cpm->freq[i] <= 0) cpm->freq[i] = 1;

			cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
		}

		vlc = NULL;
	}

	for (i = 0; i < dec->num_class; i++) {
		dec->mtfbuf[i] = i;
	}

	for (y = 0; y < dec->height; y += blksize) {
		for (x = 0; x < dec->width; x += blksize) {
			decode_qtindex(fp, dec, vlc, cpm, y, x, blksize, dec->width, level);
		}
	}

	if (dec->f_huffman == 1) {
		free_vlc(vlc);
	}
	else {
		free(cpm->freq);
	}
}

int calc_udec(DECODER *dec, int y, int x) {
	int u, k, *err_p;
	int *intra_wt_p, *back_wt_p = NULL, *for_wt_p = NULL;
	int *mi_up_wt_p = NULL, *mi_left_wt_p = NULL, *mi_ldiag_wt_p = NULL, *mi_rdiag_wt_p = NULL;
	int *intra_roff_p = NULL, *back_roff_p = NULL, *for_roff_p = NULL;
	int *mi_up_roff_p = NULL, *mi_left_roff_p = NULL, *mi_ldiag_roff_p = NULL, *mi_rdiag_roff_p = NULL;

	err_p = &dec->err[1][y][x];

	intra_roff_p = dec->intra_roff[y][x];
	intra_wt_p = dec->intra_ctx_weight;

	if (dec->back_prd_order > 0) {
		back_roff_p = dec->back_roff[y][x];
		back_wt_p = dec->back_ctx_weight;
	}

	if (dec->for_prd_order > 0) {
		for_roff_p = dec->for_roff[y][x];
		for_wt_p = dec->for_ctx_weight;
	}

	if (dec->mi_prd_order[UP] > 0) {
		mi_up_roff_p = dec->mi_up_roff[y][x];

		mi_up_wt_p = (y < dec->mi_size) ? dec->mi_borders_ctx_weight : dec->mi_up_ctx_weight;
	}

	if (dec->mi_prd_order[LEFT] > 0) {
		mi_left_roff_p = dec->mi_left_roff[y][x];

		mi_left_wt_p = (x < dec->mi_size) ? dec->mi_borders_ctx_weight : dec->mi_left_ctx_weight;
	}

	if (dec->mi_prd_order[LDIAG] > 0) {
		mi_ldiag_roff_p = dec->mi_ldiag_roff[y][x];

		mi_ldiag_wt_p = (y < dec->mi_size || x < dec->mi_size) ? dec->mi_borders_ctx_weight : dec->mi_ldiag_ctx_weight;
	}

	if (dec->mi_prd_order[RDIAG] > 0) {
		mi_rdiag_roff_p = dec->mi_rdiag_roff[y][x];
		mi_rdiag_wt_p = (y < dec->mi_size || x > dec->width - dec->mi_size) ? dec->mi_borders_ctx_weight : dec->mi_rdiag_ctx_weight;
	}

	u = 0;

	if (y == 0 && x == 0) {
		for (k = 0; k < dec->prd_order; k++) {
			u += ((dec->maxval + 1) >> 2) * (*intra_wt_p++);
		}
	}
	else {
		for (k = 0; k < dec->prd_order; k++) {
			u += err_p[*intra_roff_p++] * (*intra_wt_p++);
		}
	}

	// + dec->width needed due to extra line in the encoder
	for (k = 0; k < dec->back_prd_order; k++) {
		u += err_p[*back_roff_p++ + dec->width] * (*back_wt_p++);
	}

	// - dec->width needed due to extra line in the encoder
	for (k = 0; k < dec->for_prd_order; k++) {
		u += err_p[*for_roff_p++ - dec->width] * (*for_wt_p++);
	}

	if (y == 0 && x == 0) {
		for (k = 0; k < dec->mi_prd_order[UP]; k++) {
			u += ((dec->maxval + 1) >> 2) * (*mi_up_wt_p++);
		}
	}
	else {
		for (k = 0; k < dec->mi_prd_order[UP]; k++) {
			u += err_p[*mi_up_roff_p++] * (*mi_up_wt_p++);
		}
	}

	if (y == 0 && x == 0) {
		for (k = 0; k < dec->mi_prd_order[LEFT]; k++) {
			u += ((dec->maxval + 1) >> 2) * (*mi_left_wt_p++);
		}
	}
	else {
		for (k = 0; k < dec->mi_prd_order[LEFT]; k++) {
			u += err_p[*mi_left_roff_p++] * (*mi_left_wt_p++);
		}
	}

	if (y == 0 && x == 0) {
		for (k = 0; k < dec->mi_prd_order[LDIAG]; k++) {
			u += ((dec->maxval + 1) >> 2) * (*mi_ldiag_wt_p++);
		}
	}
	else {
		for (k = 0; k < dec->mi_prd_order[LDIAG]; k++) {
			u += err_p[*mi_ldiag_roff_p++] * (*mi_ldiag_wt_p++);
		}
	}

	if (y == 0 && x == 0) {
		for (k = 0; k < dec->mi_prd_order[RDIAG]; k++) {
			u += ((dec->maxval + 1) >> 2) * (*mi_rdiag_wt_p++);
		}
	}
	else {
		for (k = 0; k < dec->mi_prd_order[RDIAG]; k++) {
			u += err_p[*mi_rdiag_roff_p++] * (*mi_rdiag_wt_p++);
		}
	}

	u >>= 6;

	if (u > MAX_UPARA) u = MAX_UPARA;

	return (u);
}

int calc_prd(img_t ***org, DECODER *dec, int cl, int y, int x) {
	int k, prd, *coef_p;
	img_t *org_p;
	int *intra_roff_p = NULL, *back_roff_p = NULL, *for_roff_p = NULL;
	int *mi_up_roff_p = NULL, *mi_left_roff_p = NULL, *mi_ldiag_roff_p = NULL, *mi_rdiag_roff_p = NULL;

	intra_roff_p = dec->intra_roff[y][x];
	if (dec->back_prd_order > 0) back_roff_p = dec->back_roff[y][x];
	if (dec->for_prd_order > 0) for_roff_p = dec->for_roff[y][x];
	if (dec->mi_prd_order[UP] > 0) mi_up_roff_p = dec->mi_up_roff[y][x];
	if (dec->mi_prd_order[LEFT] > 0) mi_left_roff_p = dec->mi_left_roff[y][x];
	if (dec->mi_prd_order[LDIAG] > 0) mi_ldiag_roff_p = dec->mi_ldiag_roff[y][x];
	if (dec->mi_prd_order[RDIAG] > 0) mi_rdiag_roff_p = dec->mi_rdiag_roff[y][x];

	org_p = &org[1][y][x];
	coef_p = dec->predictor[cl];

	prd = 0;

	if (y == 0 && x == 0) {
		for (k = 0; k < dec->prd_order; k++) {
			prd += ((dec->maxval + 1) >> 1) * (*coef_p++);
		}
	}
	else {
		for (k = 0; k < dec->prd_order; k++) {
			prd += org_p[*intra_roff_p++] * (*coef_p++);
		}
	}

	// + dec->width needed due to extra line in the encoder
	for (k = 0; k < dec->back_prd_order; k++) {
		prd += org_p[*back_roff_p++ + dec->width] * (*coef_p++);
	}

	// - dec->width needed due to extra line in the encoder
	for (k = 0; k < dec->for_prd_order; k++) {
		prd += org_p[*for_roff_p++ - dec->width] * (*coef_p++);
	}

	if (y == 0 && x == 0) {
		for (k = 0; k < dec->mi_prd_order[UP]; k++) {
			prd += ((dec->maxval + 1) >> 1) * (*coef_p++);
		}
	}
	else {
		for (k = 0; k < dec->mi_prd_order[UP]; k++) {
			prd += org_p[*mi_up_roff_p++] * (*coef_p++);
		}
	}

	if (y == 0 && x == 0) {
		for (k = 0; k < dec->mi_prd_order[LEFT]; k++) {
			prd += ((dec->maxval + 1) >> 1) * (*coef_p++);
		}
	}
	else {
		for (k = 0; k < dec->mi_prd_order[LEFT]; k++) {
			prd += org_p[*mi_left_roff_p++] * (*coef_p++);
		}
	}

	if (y == 0 && x == 0) {
		for (k = 0; k < dec->mi_prd_order[LDIAG]; k++) {
			prd += ((dec->maxval + 1) >> 1) * (*coef_p++);
		}
	}
	else {
		for (k = 0; k < dec->mi_prd_order[LDIAG]; k++) {
			prd += org_p[*mi_ldiag_roff_p++] * (*coef_p++);
		}
	}

	if (y == 0 && x == 0) {
		for (k = 0; k < dec->mi_prd_order[RDIAG]; k++) {
			prd += ((dec->maxval + 1) >> 1) * (*coef_p++);
		}
	}
	else {
		for (k = 0; k < dec->mi_prd_order[RDIAG]; k++) {
			prd += org_p[*mi_rdiag_roff_p++] * (*coef_p++);
		}
	}

	if (prd < 0) prd = 0;
	else if (prd > dec->maxprd) prd = dec->maxprd;

	return (prd);
}

IMAGE *decode_image(FILE *fp, IMAGE *video[3], DECODER *dec) {
	int x, y, cl, gr, prd, u, e, E, p;

	// TODO: consider add org to decoder
	img_t ***org = (img_t ***) alloc_3d_array(dec->height, dec->width, 3, sizeof(img_t));

	for (y = 0; y < dec->height; y++) {
		for (x = 0; x < dec->width; x++) {
			if (dec->back_prd_order > 0) {
				org[0][y][x] = video[0]->val[y][x];
			}

			if (dec->for_prd_order > 0) {
				org[2][y][x] = video[2]->val[y][x];
			}
		}
	}

	if (dec->f_huffman == 1) {
		VLC *vlc;
		dec->vlcs = init_vlcs(dec->pmodels, dec->num_group, 1);

		for (y = 0; y < dec->height; y++) {
			for (x = 0; x < dec->width; x++) {
				cl = dec->class[y][x];
				u = calc_udec(dec, y, x);
				gr = dec->uquant[cl][u];
				prd = calc_prd(org, dec, cl, y, x);
				prd >>= (dec->coef_precision - 1);
				p = (prd + 1) >> 1;
				vlc = &dec->vlcs[gr][0];
				dec->err[1][y][x] = E = decode_vlc(fp, vlc);
				e = E2e(E, p, prd & 1, dec->maxval);
				org[1][y][x] = (img_t) (p + e);
			}
		}
	}
	else {
		PMODEL *pm;
		if (dec->pm_accuracy < 0) {
			for (y = 0; y < dec->height; y++) {
				for (x = 0; x < dec->width; x++) {
					cl = dec->class[y][x];
					u = calc_udec(dec, y, x);
					gr = dec->uquant[cl][u];
					prd = calc_prd(org, dec, cl, y, x);
					prd >>= (dec->coef_precision - 1);
					p = (prd + 1) >> 1;
					pm = dec->pmodels[gr][0];
					dec->err[1][y][x] = E = rc_decode(fp, dec->rc, pm, 0, pm->size);
					e = E2e(E, p, prd & 1, dec->maxval);
					org[1][y][x] = (img_t) (p + e);
				}
			}
		}
		else {
			int mask, shift, base;
			mask = (1 << dec->pm_accuracy) - 1;
			shift = dec->coef_precision - dec->pm_accuracy;

			for (y = 0; y < dec->height; y++) {
				for (x = 0; x < dec->width; x++) {
					cl = dec->class[y][x];
					u = calc_udec(dec, y, x);
					gr = dec->uquant[cl][u];
					prd = calc_prd(org, dec, cl, y, x);
					base = (dec->maxprd - prd + (1 << shift) / 2) >> shift;
					pm = dec->pmodels[gr][0] + (base & mask);
					base >>= dec->pm_accuracy;
					p = rc_decode(fp, dec->rc, pm, base, base+dec->maxval+1) - base;
					org[1][y][x] = (img_t) p;
					prd >>= (dec->coef_precision - 1);
					e = (p << 1) - prd;
					dec->err[1][y][x] = (e > 0)? (e - 1) : (-e);
				}
			}
		}
	}

	video[1] = alloc_image(dec->width, dec->height, dec->maxval);

	for (y = 0; y < dec->height; y++) {
		for (x = 0; x < dec->width; x++) {
			video[1]->val[y][x] = org[1][y][x];
		}
	}

	free(org);

	return (video[1]);
}

/*------------------------- decode_lookuptable -------------------------*
 |  Function decode_lookuptable
 |
 |  Purpose: Recovers lookup table from file
 |
 |  Parameters:
 |		fp				--> File to read from (IN)
 |		hist_bytes		--> Number of bytes to read from file (IN)
 |		depth			--> Length of the histogram (IN)
 |
 |  Returns:  IMGAE		--> Returns the reconstructed image
 *----------------------------------------------------------------------*/
int *decode_lookuptable(FILE *fp, int hist_bytes, int depth) {
	char *new_str = (char *) alloc_mem(sizeof(char));
	new_str[0] = '\0';

	char phrase[25] = {};
	int i = 0, j = 0, size, count = 0;

	int *forward_table = (int *) alloc_mem(sizeof(int) * (depth + 1));
	for (i = 0; i < depth; i++) {
		forward_table[i] = 0;
	}

	// Read values from file
	for (i = 0; i < hist_bytes; i++) {
		new_str = cat_str(new_str, int2bin(getbits(fp, 8), 8), 1);
	}

	size = (int) strlen(new_str);

	for (i = 0; i < size; i++) {
		phrase[0] = '\0';

		for (j = 0; j < 7; j++) {
			phrase[j] = new_str[j + i + 1];
		}

		phrase[j] = '\0';

		if (new_str[i] == '0') {
			if (strtol(phrase, NULL, 2) < 126) {
				count += strtol(phrase, NULL, 2) + 1;

				forward_table[count] = 1;

				count++;

				i = i + 7;
			}
			else {
				phrase[0] = '\0';

				if (new_str[i + 7] == '0') {
					for (j = 0; j < 8; j++) {
						phrase[j] = new_str[j + i + 8];
					}

					phrase[j] = '\0';

					count += strtol(phrase, NULL, 2) + 127;

					forward_table[count] = 1;

					count++;

					i = i + 15;
				}
				else {
					for (j = 0; j < 16; j++) {
						phrase[j] = new_str[j + i + 8];
					}

					phrase[j] = '\0';

					count += strtol(phrase, NULL, 2) + 383;

					forward_table[count] = 1;

					count++;

					i = i + 23;
				}
			}
		}
		else {
			if (strtol(phrase, NULL, 2) < 126) {
				for (j = count; j < count + strtol(phrase, NULL, 2) + 1; j++) {
					forward_table[j] = 1;
				}

				count += strtol(phrase, NULL, 2) + 2;

				i = i + 7;
			}
			else {
				phrase[0] = '\0';

				if (new_str[i + 7] == '0') {
					for (j = 0; j < 8; j++) {
						phrase[j] = new_str[j + i + 8];
					}

					phrase[j] = '\0';

					for (j = count; j < count + strtol(phrase, NULL, 2) + 127; j++) {
						forward_table[j] = 1;
					}

					count += strtol(phrase, NULL, 2) + 128;

					i = i + 15;
				}
				else {
					for (j = 0; j < 16; j++) {
						phrase[j] = new_str[j + i + 8];
					}

					phrase[j] = '\0';

					for (j = count; j < count + strtol(phrase, NULL, 2) + 383; j++) {
						forward_table[j] = 1;
					}

					count += strtol(phrase, NULL, 2) + 384;

					i = i + 23;
				}
			}
		}
	}

	free(new_str);

	count = 0;
	j = 0;

	for (i = 0; i < depth; i++) {
		count += forward_table[i];
	}

	int *backward_table = (int *) alloc_mem(sizeof(int) * count);

	for (i = 0; i < depth; i++) {
		if (forward_table[i] != 0) {
			backward_table[j] = i;
			j++;
		}
	}

	safefree((void **) &forward_table);

	return backward_table;
}

/*------------------------- histogram_unpacking -------------------------*
 |  Function histogram_unpacking
 |
 |  Purpose: Peforms the histogram packing and checks if the
 |			 total variation is lower
 |
 |  Parameters:
 |		img				--> Image to pack (IN)
 |		backward_table	--> Lookup table to use (IN)
 |
 |  Returns:  IMGAE		--> Returns the reconstructed image
 *----------------------------------------------------------------------*/
IMAGE *histogram_unpacking(IMAGE *img, const int *backward_table) {
	int x, y;
	IMAGE *u_img = alloc_image(img->width, img->height, img->maxval);

	// Perform the histogram packing of the image
	for (y = 0; y < img->height; y++) {
		for (x = 0; x < img->width; x++) {
			u_img->val[y][x] = (img_t) backward_table[img->val[y][x]];
		}
	}

	return u_img;
}

int main(int argc, char **argv) {
	int i, f, **error = NULL;
	int width, height, maxval, frames, depth, bframes, bidirectional = 0, num_comp, num_group, endianness = LITTLE_ENDIANNESS;
	int num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta, hevc, hist_bytes;
	int *backward_table = NULL;
	int prd_order[6] = {0, 0, 0, 0, 0, 0};
	int mi_prd_order[4] = {0, 0, 0, 0};
	int mi_size = 0;

	IMAGE *video[3] = {NULL, NULL, NULL};
	DECODER *dec = NULL;
	char *infile, *outfile;
	FILE *fp;

	cpu_time();
	setbuf(stdout, 0);
	infile = outfile = NULL;

	for (i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch(argv[i][1]) {
			case 'E':
				endianness = (int) strtol(argv[++i], NULL, 10);

				if (endianness != LITTLE_ENDIANNESS && endianness != BIG_ENDIANNESS) {
					endianness = LITTLE_ENDIANNESS;
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

	if (infile == NULL || outfile == NULL) {
		printf(BANNER"\n", MRP_VERSION);
		printf("usage: decmrp [options] infile outfile\n");
		printf("-E num		Endianness: little-endian = 0, big-endian = 1. Default: %s\n", "little-endian");
		printf("infile:     Input file\n");
		printf("outfile:    Output file\n");
		exit(0);
	}

	if (access(outfile, F_OK) != -1 ) {
		if (remove(outfile) != 0) {
			printf("Error deleting file %s.\n", outfile);
			return -1;
		}
	}

	fp = fileopen(infile, "rb");
	read_header(fp, &width, &height, &maxval, &frames, &depth, &bframes, &num_comp, &num_group, prd_order,
				&mi_size, mi_prd_order, &num_pmodel, &coef_precision, &pm_accuracy, &f_huffman, &quadtree_depth, &delta,
				&hevc, &hist_bytes);

	if (hist_bytes != 0) {
		backward_table = decode_lookuptable(fp, hist_bytes, (int) pow(2, depth));
	}

	printf("\nMRP-Video Decoder\n\n");
	// Print file characteristics to screen
	printf("%s (%dx%dx%dx%d) -> %s\n", infile, width, height, frames, depth, outfile);
	// Print coding parameters to screen
	printf("P = %d, V = %d, A = %d, D = %d\n\n", coef_precision, num_pmodel, pm_accuracy, delta);
	if (backward_table != NULL) {
		printf("Histogram packing on\n\n");
	}
	// Print prediction parameters to screen
	if (mi_size > 0) {
		if (frames == 1) {
			printf("Prediction order:\n\tFrame I: %d + %d %d %d %d\n\n", prd_order[0], mi_prd_order[UP],
				   mi_prd_order[LEFT], mi_prd_order[LDIAG], mi_prd_order[RDIAG]);
		}
		else if (bframes == 0) {
			printf("Prediction order:\n\tFrame I: %d + %d %d %d %d\n\tFrame P: %d %d\n\n", prd_order[0],
				   mi_prd_order[UP], mi_prd_order[LEFT], mi_prd_order[LDIAG], mi_prd_order[RDIAG], prd_order[1],
				   prd_order[2]);

			for (i = 3; i < 6; i++) {
				prd_order[i] = 0;
			}
		}
		else {
			printf("Number of B frames: %d\nPrediction order:\n\tFrame I: %d + %d %d %d %d\n\tFrame P: %d %d\n\tFrame B: %d %d %d\n\n",
				   bframes == 0 ? 0 : bframes - 1, prd_order[0], mi_prd_order[UP], mi_prd_order[LEFT],
				   mi_prd_order[LDIAG], mi_prd_order[RDIAG], prd_order[1], prd_order[2], prd_order[3], prd_order[4],
				   prd_order[5]);
		}
	}
	else {
		if (frames == 1) {
			printf("Prediction order:\n\tFrame I: %d\n\n", prd_order[0]);
		}
		else if (bframes == 0) {
			printf("Prediction order:\n\tFrame I: %d\n\tFrame P: %d %d\n\n", prd_order[0], prd_order[1], prd_order[2]);

			for (i = 3; i < 6; i++) {
				prd_order[i] = 0;
			}
		}
		else {
			printf("Number of B frames: %d\nPrediction order:\n\tFrame I: %d\n\tFrame P: %d %d\n\tFrame B: %d %d %d\n\n",
				   bframes == 0 ? 0 : bframes - 1, prd_order[0], prd_order[1], prd_order[2], prd_order[3], prd_order[4],
				   prd_order[5]);
		}
	}

	int back_reference = 0, for_reference = 0;
	int **back_ref_error = NULL, **for_ref_error = NULL;
	IMAGE **seq = (IMAGE **) alloc_mem(frames * sizeof(IMAGE));
	int ***keep_error = NULL;
	int first_frame = 0, conta = 0, final = 0;
	int **bref = NULL;

	if (hevc == 0) {
		back_reference = 0;
		for_reference = bframes;
	}
	else if (hevc == 1) {
		keep_error = (int ***)alloc_mem((bframes + 1) * sizeof(int **));

		for (f = 0; f < bframes + 1; f++) keep_error[f] = NULL;

		bref = select_bref(bframes);

		final = ((frames - 1) % bframes != 0) ? bframes * ((frames - 1) / bframes) + 1 : frames;
	}

	//Loop to decode all the frames
	f = 0;

	while(f < frames) {
		printf("Decoding frame: %03d", f);

		// Creates new DECODER structure
		if (f == 0) {
			dec = init_decoder(fp, NULL, NULL, width, height, maxval, num_comp, num_group, prd_order[0], 0, 0, mi_size,
							   mi_prd_order, num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta,
							   depth);
		}
		else if (bidirectional == 0) {
			video[0] = video[1];

			dec = init_decoder(fp, error, NULL, width, height, maxval, num_comp, num_group, prd_order[1], prd_order[2],
							   0, mi_size, mi_prd_order, num_pmodel, coef_precision, pm_accuracy, f_huffman,
							   quadtree_depth, delta, depth);

			free(error);
		}
		else {
			if ((f - first_frame) % bframes == 0) {
				dec = init_decoder(fp, back_ref_error, NULL, width, height, maxval, num_comp, num_group, prd_order[1],
								   prd_order[2], 0, mi_size, mi_prd_order, num_pmodel, coef_precision, pm_accuracy,
								   f_huffman, quadtree_depth, delta, depth);
			}
			else {
				dec = init_decoder(fp, back_ref_error, for_ref_error, width, height, maxval, num_comp, num_group,
								   prd_order[3], prd_order[4], prd_order[5], mi_size, mi_prd_order, num_pmodel,
								   coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta, depth);
			}
		}

		decode_class(fp, dec);
		decode_predictor(fp, dec);
		decode_threshold(fp, dec);

		dec->pmodels = init_pmodels(dec->num_group, dec->num_pmodel, dec->pm_accuracy, dec->pm_idx, dec->sigma, dec->maxval + 1);

		video[1] = decode_image(fp, video, dec);

		seq[f] = copy_yuv(video[1]);

		printf(" --> Process completed\n");

		if (bidirectional == 0) {
			error = get_dec_err(dec, 1);

			if (f > 0) {
				free(video[0]->val);
				free(video[0]);
			}

			f++;
		}
		else {
			if (hevc == 0) {
				if (f == 0) {
					f = bframes;

					back_ref_error = get_dec_err(dec, 1);
					video[0] = video[1];
				}
				else if (f == for_reference) {
					if (f == back_reference + 1) {
						f = frames;
					}
					else {
						f = back_reference + 1;

						for_ref_error = get_dec_err(dec, 1);

						video[2] = video[1];
					}
				}
				else if (back_reference < f && f < for_reference - 1) {
					f++;

					safefree_yuv(&video[1]);
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

					if (video[1] != video[0] && video[1] != video[2]) {
						safefree_yuv(&video[1]);
					}
					safefree_yuv(&video[0]);
					video[0] = video[2];
				}
			}
			else if (hevc == 1) {
				safefree_yuv(&video[1]);

				if (f == 0) {
					f = bref[conta][0];

					keep_error[0] = get_dec_err(dec, 1);
					for_ref_error = NULL;
					video[2] = NULL;
					back_ref_error = keep_error[0];

					video[0] = seq[f + bref[conta][1]];

					conta++;
				}
				else {
					if (conta + 1 <= bframes) {
						keep_error[f - first_frame] = get_dec_err(dec, 1);

						if (f + bref[conta][0] > 2 * first_frame && first_frame != 0) {
							f = first_frame + bref[conta][0];
						}
						else {
							f = f + bref[conta][0];
						}

						video[0] = seq[f + bref[conta][1]];

						if (conta < bframes + 1) back_ref_error = keep_error[bref[conta][3]];

						if (bref[conta][2] != -1) {
							video[2] = seq[f + bref[conta][2]];

							if (conta < bframes + 1) for_ref_error = keep_error[bref[conta][4]];
						}
						else {
							for_ref_error = NULL;
						}

						conta++;
					}
					else if (conta + 1 > bframes) {
						if (first_frame + 2 * bframes >= frames) {
							if (f + 1 < final && final != frames) {
								int numbs = frames - final;
								int aux_bframes = bframes;
								bframes = numbs;

								conta = 1;

								video[0] = seq[f + 1];

								first_frame = f + 1;
								f = frames - 1;
								final = frames;

								if (numbs == 1) {
									conta = bframes;

									video[0] = seq[frames - 2];

									first_frame = frames - 2;
								}
								else {
									bref = select_bref(bframes);
								}

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

							video[0] = seq[f + bref[0][1]];

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
		}

		free_decoder(dec);
		dec = NULL;
	}

	if (dec != NULL) free_decoder(dec);

	for (f = 0; f < frames; f++) {
		if (backward_table != NULL) {
			IMAGE *unpacked = histogram_unpacking(seq[f], backward_table);

			write_yuv(unpacked, outfile, depth, endianness);

			safefree_yuv(&unpacked);
		}
		else {
			write_yuv(seq[f], outfile, depth, endianness);
		}

		free(seq[f]->val);
		free(seq[f]);
	}

	free(seq);

	if (bidirectional == 0) {
		free(video[1]->val);
		free(video[1]);

		free(error);
	}
	else {
		if (hevc == 0) {
			if (video[1] != video[0] && video[1] != video[2]) {
				if (video[1] != NULL) safefree_yuv(&video[1]);
			}
			if (video[0] != video[2]) {
				safefree_yuv(&video[0]);
			}
			safefree_yuv(&video[2]);

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

		free(bref);
	}

	safefree((void **) &backward_table);

	fclose(fp);

	printf("\ncpu time: %.2f sec.\n", cpu_time());

	return(0);
}
