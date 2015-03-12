#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "mrp.h"

extern POINT dyx[];
extern POINT idyx[];
extern double sigma_h[], sigma_a[];
extern double qtree_prob[];
extern int bref1[][5], bref2[][5], bref3[][5], bref4[][5], bref5[][5], bref6[][5], bref7[][5], bref8[][5], bref9[][5];

/* Alternative version for 'free()' */
void safefree(void **pp){
    /* in debug mode, abort if pp is NULL */
    assert(pp);
    if (pp != NULL) {               /* safety check */
        free(*pp);                  /* deallocate chunk, note that free(NULL) is valid */
        *pp = NULL;                 /* reset original pointer */
    }
}

void safefree_yuv(IMAGE **pp){
    /* in debug mode, abort if pp is NULL */
    assert(pp);
    if (pp != NULL) {               /* safety check */
    	free((*pp)->val);
        free(*pp);                  /* deallocate chunk, note that free(NULL) is valid */
        *pp = NULL;                 /* reset original pointer */
    }
}

uint getbits(FILE *fp, int n){
	static int bitpos = 0;
	static uint bitbuf = 0;
	int x = 0;

	if (n <= 0) return (0);

	while (n > bitpos){
		n -= bitpos;
		x = (x << bitpos) | bitbuf;
		bitbuf = getc(fp) & 0xff;
		bitpos = 8;
	}

	bitpos -= n;
	x = (x << n) | (bitbuf >> bitpos);
	bitbuf &= ((1 << bitpos) - 1);

	return (x);
}

int read_class(FILE *fp){
	return(getbits(fp, 8));
}

DECODER *init_decoder(FILE *fp, int **back_ref_error, int **for_ref_error, int version, int width, int height, int maxval, int num_comp, int num_group, int prd_order, int back_prd_order, int for_prd_order, int num_pmodel, int coef_precision, int pm_accuracy, int f_huffman, int quadtree_depth, int delta){
	DECODER *dec;
	int i;

	dec = (DECODER *)alloc_mem(sizeof(DECODER));

	dec->version = version;
	dec->width = width;
	dec->height = height;
	dec->maxval = maxval;
	dec->num_comp = num_comp;
	dec->num_group = num_group;
	dec->prd_order = prd_order;
	dec->back_prd_order = back_prd_order;
	dec->for_prd_order = for_prd_order;
	dec->num_pmodel = num_pmodel;
	dec->coef_precision = coef_precision;
	dec->pm_accuracy = pm_accuracy;
	dec->f_huffman = f_huffman;
	dec->quadtree_depth = quadtree_depth;
	dec->maxprd = dec->maxval << dec->coef_precision;
	dec->delta = delta;

	dec->num_class = read_class(fp);

	dec->predictor = (int **)alloc_2d_array(dec->num_class, dec->prd_order + dec->back_prd_order + dec->for_prd_order, sizeof(int));

	dec->err = (int ***)alloc_3d_array(dec->height, dec->width, 3, sizeof(int));

	if(back_ref_error != NULL){
		int y, x;

		for(y = 0; y < dec->height; y++){
			for(x = 0; x < dec->width; x++){
				dec->err[0][y][x] = back_ref_error[y][x];
			}
		}
	}

	if(for_ref_error != NULL){
		int y, x;

		for(y = 0; y < dec->height; y++){
			for(x = 0; x < dec->width; x++){
				dec->err[2][y][x] = for_ref_error[y][x];
			}
		}
	}

	dec->ctx_weight = init_ctx_weight(dec->prd_order, dec->back_prd_order, dec->for_prd_order, dec->delta);

	// Quadtree map
	if (dec->quadtree_depth > 0){
		int x, y, xx, yy;

		yy = (dec->height + MAX_BSIZE - 1) / MAX_BSIZE;
		xx = (dec->width + MAX_BSIZE - 1) / MAX_BSIZE;

		for (i = dec->quadtree_depth - 1; i >= 0; i--){
			dec->qtmap[i] = (char **)alloc_2d_array(yy, xx, sizeof(char));

			for (y = 0; y < yy; y++){
				for (x = 0; x < xx; x++){
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

	if (dec->num_pmodel > 1){
		dec->pm_idx = (int *)alloc_mem(dec->num_group * sizeof(int));
	}
	else{
		dec->pm_idx = NULL;
	}

	dec->spm.freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
	dec->spm.cumfreq = &(dec->spm.freq[MAX_SYMBOL]);

	//Huffman coding
	if (dec->f_huffman == 1){
		dec->sigma = sigma_h;
	}
	else{
		dec->sigma = sigma_a;
	}

	dec->mtfbuf = (int *)alloc_mem(dec->num_class * sizeof(int));

	return (dec);
}

int **get_dec_err(DECODER *dec, int pos){
	int y, x;
	int **error = (int **)alloc_2d_array(dec->height + 1, dec->width, sizeof(int));

	for (y = 0; y < dec->height; y++){
		for (x = 0; x < dec->width; x++){
			error[y][x] = dec->err[pos][y][x];
		}
	}

	return error;
}


void free_decoder(DECODER *dec){
	int i, j, gr, num_subpm;

	free(dec->predictor);
	free(dec->err);
	free(dec->ctx_weight);

	if (dec->quadtree_depth > 0){
		for (i = dec->quadtree_depth - 1; i >= 0; i--){
			free(dec->qtmap[i]);
		}
	}

	free(dec->class);
	free(dec->uquant);

	if (dec->num_pmodel > 1){
		free(dec->pm_idx);
	}

	free(dec->spm.freq);

	if(dec->pm_accuracy < 0){
		num_subpm = 1;
	}
	else{
		num_subpm = 1 << dec->pm_accuracy;
	}
	for (gr = 0; gr < dec->num_group; gr++){
		for (j = 0; j < num_subpm; j++){
			PMODEL *aux = &dec->pmodels[gr][0][j];
			free(aux->freq);
		}
	}
	free(dec->pmodels[0][0]);
	free(dec->pmodels);

	free(dec->mtfbuf);

	if (dec->f_huffman == 0){
		free(dec->rc);
	}

	free(dec);
}

void read_header(FILE *fp, int *version, int *width, int *height, int *maxval, int *frames, int *bframes, int *num_comp, int *num_group, int *prd_order, int *num_pmodel, int *coef_precision, int *pm_accuracy, int *f_huffman, int *quadtree_depth, int *delta, int *diff, int *hevc){
	int i = 0;

	if (getbits(fp, 16) != MAGIC_NUMBER){
		fprintf(stderr, "Not a compressed file!\n");
		exit(1);
	}

	*version = getbits(fp, 8);
	*width = getbits(fp, 16);
	*height = getbits(fp, 16);
	*maxval = getbits(fp, 16);
	*frames = getbits(fp, 16);
	*bframes = getbits(fp, 8);
	*hevc = getbits(fp, 1);
	*num_comp = getbits(fp, 4);
	*num_group = getbits(fp, 6);
	for(i = 0; i < 6; i++){
		prd_order[i] = getbits(fp, 8);
	}
	*diff = getbits(fp, 1);
	*delta = getbits(fp, 8);
	*num_pmodel = getbits(fp, 6) + 1;
	*coef_precision = getbits(fp, 4) + 1;
	*pm_accuracy = getbits(fp, 3) - 1;
	*f_huffman = getbits(fp, 1);
	*quadtree_depth = (getbits(fp, 1))? QUADTREE_DEPTH : -1;
	getbits(fp, 5);
}

int decode_vlc(FILE *fp, VLC *vlc){
	int i, k, min, off;
	uint code;

	code = min = off = k = 0;

	for (i = 0; i < vlc->max_len; i++){
		code = (code << 1) | getbits(fp, 1);
		k = vlc->off[i];

		if (k < 0){
			min <<= 1;
		}
		else{
			if (code <= vlc->code[vlc->index[k]]) break;

			min = (vlc->code[vlc->index[k]] + 1) << 1;
			off = k + 1;
		}
	}

	i = off + code - min;

	return (vlc->index[i]);
}

int decode_golomb(FILE *fp, int m){
	int v = 0;

	while (getbits(fp, 1) == 0){
		v++;
	}

	v = (v << m) | getbits(fp, m);

	return (v);
}

void decode_predictor(FILE *fp, DECODER *dec){
	int k, m, cl, coef, sgn;

	int prd_order = dec->prd_order + dec->back_prd_order + dec->for_prd_order;

	if (dec->f_huffman == 1){
		for (k = 0; k < prd_order; k++){
			m = getbits(fp, 4);

			for (cl = 0; cl < dec->num_class; cl++){
				coef = decode_golomb(fp, m);

				if (coef > 0){
					sgn = getbits(fp, 1);

					if (sgn){
						coef = -coef;
					}
				}

				dec->predictor[cl][k] = coef;
			}
		}
	}
	else{
		PMODEL *pm;

		pm = &dec->spm;
		pm->size = MAX_COEF + 18;
		pm->cumfreq[MAX_COEF + 2] = 0;

		for(k = MAX_COEF + 2; k < pm->size; k++){
			pm->freq[k] = 1;
			pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
		}

		for (k = 0; k < prd_order; k++){
			m = rc_decode(fp, dec->rc, pm, MAX_COEF + 2, MAX_COEF + 18) - (MAX_COEF + 2);
			set_spmodel(pm, MAX_COEF + 1, m);

			for (cl = 0; cl < dec->num_class; cl++){
				coef = rc_decode(fp, dec->rc, pm, 0, MAX_COEF + 1);

				if (coef > 0){
					sgn = rc_decode(fp, dec->rc, pm, MAX_COEF+2, MAX_COEF+4) - (MAX_COEF + 2);

					if (sgn){
						coef = -coef;
					}
				}

				dec->predictor[cl][k] = coef;
			}
		}
	}

	return;
}

void decode_threshold(FILE *fp, DECODER *dec){
	int cl, gr, m, u, k;

	if (dec->f_huffman == 1){
		m = getbits(fp, 4);

		for (cl = 0; cl < dec->num_class; cl++){
			k = u = 0;

			for (gr = 0; gr < dec->num_group; gr++){
				if (k > MAX_UPARA || gr == dec->num_group - 1){
					k = MAX_UPARA + 1;
				}
				else{
					if (getbits(fp, 1)) k += decode_golomb(fp, m) + 1;
				}

				for (; u < k; u++) dec->uquant[cl][u] = gr;
			}
		}

		if (dec->num_pmodel > 1){
			for (k = 1; (1 << k) < dec->num_pmodel; k++);

			for (gr = 0; gr < dec->num_group; gr++){
				dec->pm_idx[gr] = getbits(fp, k);
			}
		}
	}
	else{
		PMODEL *pm;

		pm = &dec->spm;
		pm->size = 16;
		pm->cumfreq[0] = 0;

		for (k = 0; k < pm->size; k++){
			pm->freq[k] = 1;
			pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
		}

		m = rc_decode(fp, dec->rc, pm, 0, pm->size);
		set_spmodel(pm, MAX_UPARA + 2, m);

		for (cl = 0; cl < dec->num_class; cl++){
			k = u = 0;

			for (gr = 0; gr < dec->num_group; gr++){
				if (k > MAX_UPARA || gr == dec->num_group - 1){
					k = MAX_UPARA + 1;
				}
				else{
					k += rc_decode(fp, dec->rc, pm, 0, pm->size - k);
				}

				for (; u < k; u++) dec->uquant[cl][u] = gr;
			}
		}

		if (dec->num_pmodel > 1){
			pm->size = dec->num_pmodel;
			pm->freq[0] = 0;

			for (k = 0; k < pm->size; k++){
				pm->freq[k] = 1;
				pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
			}

			for (gr = 0; gr < dec->num_group; gr++){
				dec->pm_idx[gr] = rc_decode(fp, dec->rc, pm, 0, pm->size);
			}
		}
	}

	return;
}

void decode_qtindex(FILE *fp, DECODER *dec, VLC *vlc, PMODEL *cpm, int tly, int tlx, int blksize, int width, int level){
	int i, cl, y, x, bry, brx, ctx;
	char **qtmap;
	PMODEL *pm;

	brx = (tlx + blksize < dec->width) ? (tlx + blksize) : dec->width;
	bry = (tly + blksize < dec->height) ? (tly + blksize) : dec->height;

	if (tlx >= brx || tly >= bry) return;

	if (level > 0){
		ctx = 0;
		qtmap = dec->qtmap[level - 1];
		y = (tly / MIN_BSIZE) >> level;
		x = (tlx / MIN_BSIZE) >> level;

		if (y > 0){
			if (qtmap[y - 1][x] == 1) ctx++;
			if (brx < width && qtmap[y - 1][x + 1] == 1) ctx++;
		}

		if (x > 0 && qtmap[y][x - 1] == 1) ctx++;

		ctx = ((level - 1) * 4 + ctx) << 1;

		if (dec->f_huffman == 1){
			i = getbits(fp, 1);
		}
		else{
			pm = &dec->spm;
			i = rc_decode(fp, dec->rc, pm, ctx, ctx + 2) - ctx;
		}

		if (i == 1){
			qtmap[y][x] = 1;
			blksize >>= 1;
			decode_qtindex(fp, dec, vlc, cpm, tly, tlx, blksize, width, level - 1);
			decode_qtindex(fp, dec, vlc, cpm, tly, tlx + blksize, blksize, width, level - 1);
			decode_qtindex(fp, dec, vlc, cpm, tly + blksize, tlx, blksize, width, level - 1);
			decode_qtindex(fp, dec, vlc, cpm, tly + blksize, tlx + blksize, blksize, brx, level - 1);

			return;
		}
	}

	if (dec->f_huffman == 1){
		i = decode_vlc(fp, vlc);
	}
	else{
		i = rc_decode(fp, dec->rc, cpm, 0, cpm->size);
	}

	mtf_classlabel(dec->class, dec->mtfbuf, tly, tlx,blksize, width, dec->num_class);

	for (cl = 0; cl < dec->num_class; cl++){
		if (dec->mtfbuf[cl] == i) break;
	}

	for (y = tly; y < bry; y++){
		for (x = tlx; x < brx; x++){
			dec->class[y][x] = cl;
		}
	}

	return;
}

void decode_class(FILE *fp, DECODER *dec){
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
				pm->freq[(ctx << 1) + 1] = p * (1 << 10);
				p = 1.0 - p;
				pm->freq[(ctx << 1)] = p * (1 << 10);
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
			cpm->freq[i] = p * (1 << 10);

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
	else{
		free(cpm->freq);
	}

	return;
}

int calc_udec(DECODER *dec, int y, int x){
	int ry, rx, u, k;
	u = 0;

	int min_dx, max_dx, min_dy;
	int bmin_dx, bmax_dx, bmin_dy, bmax_dy;
	int fmin_dx, fmax_dx, fmin_dy, fmax_dy;

	min_dx = max_dx = min_dy = 0;
	bmin_dx = bmax_dx = bmin_dy = bmax_dy = 0;
	fmin_dx = fmax_dx = fmin_dy = fmax_dy = 0;

	//Values to check for special cases
	for (k = 0; k < dec->prd_order; k++){
		ry = dyx[k].y;
		rx = dyx[k].x;

		if (ry < min_dy) min_dy = ry;
		if (rx < min_dx) min_dx = rx;
		if (rx > max_dx) max_dx = rx;
	}

	for (k = 0; k < dec->back_prd_order - 1; k++){
		ry = idyx[k].y;
		rx = idyx[k].x;

		if (ry < bmin_dy) bmin_dy = ry;
		if (ry > bmax_dy) bmax_dy = ry;
		if (rx < bmin_dx) bmin_dx = rx;
		if (rx > bmax_dx) bmax_dx = rx;
	}

	for (k = 0; k < dec->for_prd_order - 1; k++){
		ry = idyx[k].y;
		rx = idyx[k].x;

		if (ry < fmin_dy) fmin_dy = ry;
		if (ry > fmax_dy) fmax_dy = ry;
		if (rx < fmin_dx) fmin_dx = rx;
		if (rx > fmax_dx) fmax_dx = rx;
	}

	min_dy = -min_dy;
	min_dx = -min_dx;
	max_dx = dec->width - max_dx;

	bmin_dy = -bmin_dy;
	bmin_dx = -bmin_dx;
	bmax_dx = dec->width - bmax_dx;
	bmax_dy = dec->height - bmax_dy;

	fmin_dy = -fmin_dy;
	fmin_dx = -fmin_dx;
	fmax_dx = dec->width - fmax_dx;
	fmax_dy = dec->height - fmax_dy;

	if (y >= min_dy && x >= min_dx && x <= max_dx){
		for (k = 0; k < dec->prd_order; k++){
			ry = y + dyx[k].y;
			rx = x + dyx[k].x;
			u += dec->err[1][ry][rx] * dec->ctx_weight[k];
		}
	}
	else if (y == 0){
		if (x == 0){
			for (k = 0; k < dec->prd_order; k++){
				u += ((dec->maxval + 1) >> 2) * dec->ctx_weight[k];
			}
		}
		else{
			ry = 0;

			for (k =0; k < dec->prd_order; k++){
				rx = x + dyx[k].x;

				if (rx < 0) rx = 0;
				else if (rx >= x) rx = x - 1;

				u += dec->err[1][ry][rx] * dec->ctx_weight[k];
			}
		}
	}
	else{
		if (x == 0){
			for (k = 0; k < dec->prd_order; k++){
				ry = y + dyx[k].y;

				if (ry < 0) ry = 0;
				else if (ry >= y) ry = y - 1;

				rx = x + dyx[k].x;

				if (rx < 0) rx = 0;

				u += dec->err[1][ry][rx] * dec->ctx_weight[k];
			}
		}
		else{
			for (k = 0; k < dec->prd_order; k++){
				ry = y + dyx[k].y;

				if (ry < 0) ry = 0;

				rx = x + dyx[k].x;

				if (rx < 0) rx = 0;
				else if (rx >= dec->width) rx = dec->width - 1;

				u += dec->err[1][ry][rx] * dec->ctx_weight[k];
			}
		}
	}

	//If inter prd order is different from zero and prd_order is less that NUM_UPELS
	if(dec->back_prd_order > 0){
		u += dec->err[0][y][x] * dec->ctx_weight[dec->prd_order];

		if (y >= bmin_dy && x >= bmin_dx && x <= bmax_dx && y < bmax_dy){
			for (k = 0; k < dec->back_prd_order - 1; k++){
				ry = y + idyx[k].y;
				rx = x + idyx[k].x;

				u += dec->err[0][ry][rx] * dec->ctx_weight[k + dec->prd_order + 1];
			}
		}
		else{
			for (k = 0; k < dec->back_prd_order - 1; k++){
				ry = y + idyx[k].y;
				rx = x + idyx[k].x;

				if(ry < 0 || rx < 0 || ry >= dec->height || rx >= dec->width){
					u += dec->err[0][y][x] * dec->ctx_weight[k + dec->prd_order + 1];
				}
				else{
					u += dec->err[0][ry][rx] * dec->ctx_weight[k + dec->prd_order + 1];
				}
			}
		}
	}

	//If inter prd order is different from zero and prd_order is less that NUM_UPELS
	if(dec->for_prd_order > 0){
		u += dec->err[2][y][x] * dec->ctx_weight[dec->prd_order + dec->back_prd_order];

		if (y >= fmin_dy && x >= fmin_dx && x <= fmax_dx && y < fmax_dy){
			for (k = 0; k < dec->for_prd_order - 1; k++){
				ry = y + idyx[k].y;
				rx = x + idyx[k].x;

				u += dec->err[2][ry][rx] * dec->ctx_weight[k + dec->prd_order + dec->back_prd_order + 1];
			}
		}
		else{
			for (k = 0; k < dec->for_prd_order - 1; k++){
				ry = y + idyx[k].y;
				rx = x + idyx[k].x;

				if(ry < 0 || rx < 0 || ry >= dec->height || rx >= dec->width){
					u += dec->err[2][y][x] * dec->ctx_weight[k + dec->prd_order + dec->back_prd_order + 1];
				}
				else{
					u += dec->err[2][ry][rx] * dec->ctx_weight[k + dec->prd_order + dec->back_prd_order + 1];
				}
			}
		}
	}

	u >>= 6;

	if (u > MAX_UPARA) u = MAX_UPARA;

	return (u);
}

int calc_prd(IMAGE *video[3], DECODER *dec, int cl, int y, int x){
	int k, prd, rx, ry;

	int dy, dx, min_dx, max_dx, min_dy;
	int bmin_dx, bmax_dx, bmin_dy, bmax_dy;
	int fmin_dx, fmax_dx, fmin_dy, fmax_dy;
	int min_abs_dx, max_abs_dx, min_abs_dy, max_abs_dy;

	min_dx = max_dx = min_dy = 0;
	bmin_dx = bmax_dx = bmin_dy = bmax_dy = 0;
	fmin_dx = fmax_dx = fmin_dy = fmax_dy = 0;

	//Values to check for special cases
	for (k = 0; k < dec->prd_order; k++){
		dy = dyx[k].y;
		dx = dyx[k].x;

		if (dy < min_dy) min_dy = dy;
		if (dx < min_dx) min_dx = dx;
		if (dx > max_dx) max_dx = dx;
	}

	for (k = 0; k < dec->back_prd_order - 1; k++){
		dy = idyx[k].y;
		dx = idyx[k].x;

		if (dy < bmin_dy) bmin_dy = dy;
		if (dy > bmax_dy) bmax_dy = dy;
		if (dx < bmin_dx) bmin_dx = dx;
		if (dx > bmax_dx) bmax_dx = dx;
	}

	for (k = 0; k < dec->for_prd_order - 1; k++){
		dy = idyx[k].y;
		dx = idyx[k].x;

		if (dy < fmin_dy) fmin_dy = dy;
		if (dy > fmax_dy) fmax_dy = dy;
		if (dx < fmin_dx) fmin_dx = dx;
		if (dx > fmax_dx) fmax_dx = dx;
	}

	min_abs_dy = (min_dy < bmin_dy && min_dy < fmin_dy ? -min_dy : bmin_dy < fmin_dy ? -bmin_dy : -fmin_dy);
	min_abs_dx = (min_dx < bmin_dx && min_dx < fmin_dx ? -min_dx : bmin_dx < fmin_dx ? -bmin_dx : -fmin_dx);

	if(dec->back_prd_order <= 1 && dec->for_prd_order <= 1){
		max_abs_dy = 0;
	}
	else{
		max_abs_dy = (fmax_dy > bmax_dy ? dec->height - fmax_dy : dec->height - bmax_dy);
	}

	max_abs_dx = (max_dx > bmax_dx && max_dx > fmax_dx ? dec->width - max_dx : bmax_dx > fmax_dx ? dec->width - bmax_dx : dec->width - fmax_dx);

	prd = 0;

	if (y >= min_abs_dy && x >= min_abs_dx && x < max_abs_dx && y < max_abs_dy){
		for (k = 0; k < dec->prd_order; k++){
			ry = y + dyx[k].y;
			rx = x + dyx[k].x;

			prd += dec->predictor[cl][k] * video[1]->val[ry][rx];
		}
	}
	else if (y == 0){
		if (x == 0){
			for (k = 0; k < dec->prd_order; k++){
				prd += dec->predictor[cl][k];
			}

			prd *= ((video[1]->maxval + 1) >> 1);
		}
		else{
			ry = 0;

			for (k = 0; k < dec->prd_order; k++){
				rx = x + dyx[k].x;

				if (rx < 0) rx = 0;
				else if (rx >= x) rx = x - 1;

				prd += dec->predictor[cl][k] * video[1]->val[ry][rx];
			}
		}
	}
	else{
		if (x == 0){
			for (k = 0; k < dec->prd_order; k++){
				ry = y + dyx[k].y;

				if (ry < 0) ry = 0;
				else if (ry >= y) ry = y - 1;

				rx = x + dyx[k].x;

				if (rx < 0) rx = 0;

				prd += dec->predictor[cl][k] * video[1]->val[ry][rx];
			}
		}
		else{
			for (k = 0; k < dec->prd_order; k++){
				ry = y + dyx[k].y;

				if (ry < 0) ry = 0;

				rx = x + dyx[k].x;

				if (rx < 0) rx = 0;
				else if (rx >= video[1]->width) rx = video[1]->width - 1;

				prd += dec->predictor[cl][k] * video[1]->val[ry][rx];
			}
		}
	}

	//Back Inter prediction calculation
	if(dec->back_prd_order == 1){
		prd += dec->predictor[cl][dec->prd_order] * video[0]->val[y][x];
	}
	if(dec->back_prd_order > 1){
		prd += dec->predictor[cl][dec->prd_order] * video[0]->val[y][x];

		if (y >= min_abs_dy && x >= min_abs_dx && x < max_abs_dx && y < max_abs_dy){
			for (k = 0; k < dec->back_prd_order - 1; k++){
				ry = y + idyx[k].y;
				rx = x + idyx[k].x;

				prd += dec->predictor[cl][k + dec->prd_order + 1] * video[0]->val[ry][rx];
			}
		}
		else{
			for (k = 0; k < dec->back_prd_order - 1; k++){
				ry = y + idyx[k].y;
				rx = x + idyx[k].x;

				if(ry < 0 || rx < 0 || ry >= dec->height || rx >= dec->width){
					prd += dec->predictor[cl][k + dec->prd_order + 1] * video[0]->val[y][x];
				}
				else{
					prd += dec->predictor[cl][k + dec->prd_order + 1] * video[0]->val[ry][rx];
				}
			}
		}
	}

	//Forward Inter prediction calculation
	if(dec->for_prd_order == 1){
		prd += dec->predictor[cl][dec->prd_order + dec->back_prd_order] * video[2]->val[y][x];
	}
	if(dec->for_prd_order > 1){
		prd += dec->predictor[cl][dec->prd_order + dec->back_prd_order] * video[2]->val[y][x];

		if (y >= min_abs_dy && x >= min_abs_dx && x < max_abs_dx && y < max_abs_dy){
			for (k = 0; k < dec->for_prd_order - 1; k++){
				ry = y + idyx[k].y;
				rx = x + idyx[k].x;

				prd += dec->predictor[cl][k + dec->prd_order + dec->back_prd_order + 1] * video[2]->val[ry][rx];
			}
		}
		else{
			for (k = 0; k < dec->for_prd_order - 1; k++){
				ry = y + idyx[k].y;
				rx = x + idyx[k].x;

				if(ry < 0 || rx < 0 || ry >= dec->height || rx >= dec->width){
					prd += dec->predictor[cl][k + dec->prd_order + dec->back_prd_order + 1] * video[2]->val[y][x];
				}
				else{
					prd += dec->predictor[cl][k + dec->prd_order + dec->back_prd_order + 1] * video[2]->val[ry][rx];
				}
			}
		}
	}

	if (prd < 0) prd = 0;
	else if (prd > dec->maxprd) prd = dec->maxprd;

	return (prd);
}
void print_contexts(int modo, uint cumfreq, uint freq, uint totfreq, int max, int min, int prd){
	if(modo == 0) system("rm /tmp/decoder_context.txt");

	FILE *fp;
	fp = fileopen("/tmp/decoder_context.txt", "a");

	if(modo != -1){
		fprintf(fp, "Frame: %d\n", modo);
		fprintf(fp, "Cum. Freq.\tFreq.\tTot. Freq.\tMax\tMin\tPrd:\n");
	}
	else{
		fprintf(fp, "%d\t\t%d\t%d\t\t%d\t%d\t%d\n", cumfreq, freq, totfreq, max, min, prd);
	}

	fclose(fp);
}
IMAGE *decode_image(FILE *fp, IMAGE *video[3], DECODER *dec){
	int x, y, cl, gr, prd, u, e, E, p;

	video[1] = alloc_image(dec->width, dec->height, dec->maxval);

	if (dec->f_huffman == 1){
		VLC *vlc;
		dec->vlcs = init_vlcs(dec->pmodels, dec->num_group, 1);

		for (y = 0; y < dec->height; y++){
			for (x = 0; x < dec->width; x++){
				cl = dec->class[y][x];
				u = calc_udec(dec, y, x);
				gr = dec->uquant[cl][u];
				prd = calc_prd(video, dec, cl, y, x);
				prd >>= (dec->coef_precision - 1);
				p = (prd + 1) >> 1;
				vlc = &dec->vlcs[gr][0];
				dec->err[1][y][x] = E = decode_vlc(fp, vlc);
				e = E2e(E, p, prd & 1, dec->maxval);
				video[1]->val[y][x] = p + e;
			}
		}
	}
	else{
		PMODEL *pm;
		if (dec->pm_accuracy < 0){
			for (y = 0; y < dec->height; y++){
				for (x = 0; x < dec->width; x++){
					cl = dec->class[y][x];
					u = calc_udec(dec, y, x);
					gr = dec->uquant[cl][u];
					prd = calc_prd(video, dec, cl, y, x);
					prd >>= (dec->coef_precision - 1);
					p = (prd + 1) >> 1;
					pm = dec->pmodels[gr][0];
					dec->err[1][y][x] = E = rc_decode(fp, dec->rc, pm, 0, pm->size);
					e = E2e(E, p, prd & 1, dec->maxval);
					video[1]->val[y][x] = p + e;
				}
			}
		}
		else{
			int mask, shift, base;
			mask = (1 << dec->pm_accuracy) - 1;
			shift = dec->coef_precision - dec->pm_accuracy;

			for (y = 0; y < dec->height; y++){
				for (x = 0; x < dec->width; x++){
					cl = dec->class[y][x];
					u = calc_udec(dec, y, x);
					gr = dec->uquant[cl][u];
					prd = calc_prd(video, dec, cl, y, x);
					base = (dec->maxprd - prd + (1 << shift) / 2) >> shift;
					pm = dec->pmodels[gr][0] + (base & mask);
					base >>= dec->pm_accuracy;

//					print_contexts(-1, 0, 0, pm->cumfreq[base+dec->maxval+1] - pm->cumfreq[base], base+dec->maxval+1, base, prd);

					p = rc_decode(fp, dec->rc, pm, base, base+dec->maxval+1) - base;
					video[1]->val[y][x] = p;
					prd >>= (dec->coef_precision - 1);
					e = (p << 1) - prd;
					dec->err[1][y][x] = (e > 0)? (e - 1) : (-e);
				}
			}
		}
	}

	return (video[1]);
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
IMAGE *read_yuv(char *filename, int height, int width, int frame){
	int i, j;
	IMAGE *img;
	FILE *fp;

	//Open file
	fp = fileopen(filename, "rb");

	// Check if image dimensions are correct (It has to be multiple of BASE_BSIZE)
	if ((width % BASE_BSIZE) || (height % BASE_BSIZE)){
		fprintf(stderr, "Image width and height must be multiples of %d!\n", BASE_BSIZE);
		exit(1);
	}

	// Image allocation
	img = alloc_image(width, height, 255);

	if(frame > 0) fseek(fp, img->height * img->width * 1.5 * frame, SEEK_SET);

	for (i = 0; i < img->height; i++){
		for (j = 0; j < img->width; j++){
			img->val[i][j] = (img_t)fgetc(fp);
		}
	}

	for (i = 0; i < img->height / 2; i++){
		for (j = 0; j < img->width / 2; j++){
			fgetc(fp);
		}
	}

	for (i = 0; i < img->height / 2; i++){
		for (j = 0; j < img->width / 2; j++){
			fgetc(fp);
		}
	}

	fclose(fp);
	return (img);
}

IMAGE* sum_diff(IMAGE* ref, IMAGE* diff, int frame){
	int x, y;
	// Image allocation
	IMAGE *cur = alloc_image(ref->width, ref->height, 255);

	for(y = 0; y < ref->height; y++){
		for(x = 0; x < ref->width; x++){
			cur->val[y][x] = diff->val[y][x] + ref->val[y][x] - 127;
		}
	}

	if(frame > 1){
		safefree_yuv(&ref);
	}

	return cur;
}

char *decode_extra_info(FILE *fp, int *num_pels){
	int i;
	char *extra_info = NULL;

	*num_pels = getbits(fp, 16);

	if(*num_pels != 0){
		extra_info = (char *) alloc_mem(*num_pels * sizeof(char));

		for(i = 0; i < *num_pels; i++){
			extra_info[i] = getbits(fp, 8);
		}
	}

	return extra_info;
}

int main(int argc, char **argv){
	int i, f, **error = NULL;
	int version, width, height, maxval, frames, bframes, num_comp, num_group, prd_order[6] = {0, 0, 0, 0, 0, 0}, num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta, diff, hevc;
	IMAGE *video[3] = {NULL, NULL, NULL};
	DECODER *dec;
	char *infile, *outfile;
	FILE *fp;
	cpu_time();
	setbuf(stdout, 0);
	infile = outfile = NULL;

	IMAGE *img_diff = NULL;
	int num_pels = 0, x, y;
	char *extra_info = NULL;

	for (i = 1; i < argc; i++){
		if (infile == NULL){
			infile = argv[i];
		}
		else{
			outfile = argv[i];
		}
	}

	if (infile == NULL || outfile == NULL){
		printf(BANNER"\n", 0.1 * VERSION);
		printf("usage: decmrp infile outfile\n");
		printf("infile:     Input file\n");
		printf("outfile:    Output file\n");
		exit(0);
	}

	if(access(outfile, F_OK) != -1 ){
		if(remove(outfile) != 0){
			printf("Error deleting file %s.\n", outfile);
			return -1;
		}
	}

	fp = fileopen(infile, "rb");
	read_header(fp, &version, &width, &height, &maxval, &frames, &bframes, &num_comp, &num_group, prd_order, &num_pmodel, &coef_precision, &pm_accuracy, &f_huffman, &quadtree_depth, &delta, &diff, &hevc);

	printf("\nMRP-Video Decoder\n\n");
	// Print file characteristics to screen
	printf("%s -> %s (%dx%dx%d)\n", infile, outfile, width, height, frames);
	// Print coding parameters to screen
	printf("P = %d, V = %d, A = %d, D = %d, p = %s\n\n", coef_precision, num_pmodel, pm_accuracy, delta, (diff == 1) ? "on": "off");
	printf("Number of B frames: %d\nPrediction order:\n\tFrame I: %d\n\tFrame P: %d %d\n\tFrame B: %d %d %d\n\n", bframes == 0 ? 0 : bframes - 1, prd_order[0], prd_order[1], prd_order[2], prd_order[3], prd_order[4], prd_order[5]);

	if(bframes == 0){
		for(f = 0; f < frames; f++){
			printf("Decoding frame: %03d", f);

			if(f == 0){
				dec = init_decoder(fp, NULL, NULL, version, width, height, maxval, num_comp, num_group, prd_order[0], 0, 0, num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta);
			}
			else{
				video[0] = video[1];

				dec = init_decoder(fp, error, NULL, version, width, height, maxval, num_comp, num_group, prd_order[1], prd_order[2], 0, num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta);
				free(error);
			}

			if (dec->f_huffman == 0){
				dec->rc = rc_init();
				rc_startdec(fp, dec->rc);
			}

			decode_class(fp, dec);
			decode_predictor(fp, dec);
			decode_threshold(fp, dec);

			dec->pmodels = init_pmodels(dec->num_group, dec->num_pmodel, dec->pm_accuracy, dec->pm_idx, dec->sigma, dec->maxval + 1);

			video[1] = decode_image(fp, video, dec);

			if(f > 0 && diff != 0){
				extra_info = decode_extra_info(fp, &num_pels);

				IMAGE *aux_image;
				if(num_pels != 0 && extra_info != NULL){
					aux_image = alloc_image(dec->width, dec->height, 255);
					int conta = 0;

					for(y = 0; y < dec->height; y++){
						for(x = 0; x < dec->width; x++){
							aux_image->val[y][x] = video[1]->val[y][x];

							if((video[1]->val[y][x] == 0 || video[1]->val[y][x] == 255) && conta < num_pels){
								aux_image->val[y][x] += extra_info[conta];
								conta++;
							}
						}
					}
				}

				if(extra_info != NULL) free(extra_info);

				if(num_pels == 0){
					img_diff = sum_diff(img_diff, video[1], f);
				}
				else{
					img_diff = sum_diff(img_diff, aux_image, f);
					free(aux_image);
				}

				write_yuv(img_diff, outfile);
			}
			else{
				write_yuv(video[1], outfile);
			}

			if(diff != 0 && f == 0) img_diff = video[1];

			error = get_dec_err(dec, 1);
			free_decoder(dec);

			if(f > 0){
				free(video[0]->val);
				free(video[0]);
			}

			printf(" --> Process completed\n");
		}

		if(diff == 1){
			free(img_diff->val);
			free(img_diff);
		}
		free(video[1]->val);
		free(video[1]);

		free(error);
	}
	else if(hevc == 0){
		int back_reference = 0, for_reference = bframes;
		int **back_ref_error = NULL, **for_ref_error = NULL;
		IMAGE **seq = alloc_mem(frames * sizeof(IMAGE));
		f = 0;

		while(f < frames){
		//for(f = 0; f < frames; f++){
			printf("Decoding frame: %03d", f);

			if(f == 0){
				dec = init_decoder(fp, NULL, NULL, version, width, height, maxval, num_comp, num_group, prd_order[0], 0, 0, num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta);
			}
			else if(f % bframes == 0){
				dec = init_decoder(fp, back_ref_error, NULL, version, width, height, maxval, num_comp, num_group, prd_order[1], prd_order[2], 0, num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta);
			}
			else{
				dec = init_decoder(fp, back_ref_error, for_ref_error, version, width, height, maxval, num_comp, num_group, prd_order[3], prd_order[4], prd_order[5], num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta);
			}

			if (dec->f_huffman == 0){
				dec->rc = rc_init();
				rc_startdec(fp, dec->rc);
			}

			decode_class(fp, dec);
			decode_predictor(fp, dec);
			decode_threshold(fp, dec);

			dec->pmodels = init_pmodels(dec->num_group, dec->num_pmodel, dec->pm_accuracy, dec->pm_idx, dec->sigma, dec->maxval + 1);

			video[1] = decode_image(fp, video, dec);

			if(f > 0 && diff != 0){
				extra_info = decode_extra_info(fp, &num_pels);

				IMAGE *aux_image;
				if(num_pels != 0 && extra_info != NULL){
					aux_image = alloc_image(dec->width, dec->height, 255);
					int conta = 0;

					for(y = 0; y < dec->height; y++){
						for(x = 0; x < dec->width; x++){
							aux_image->val[y][x] = video[1]->val[y][x];

							if((video[1]->val[y][x] == 0 || video[1]->val[y][x] == 255) && conta < num_pels){
								aux_image->val[y][x] += extra_info[conta];
								conta++;
							}
						}
					}
				}

				if(extra_info != NULL) free(extra_info);

				if(num_pels == 0){
					seq[f] = copy_yuv(video[1]);
				}
				else{
					seq[f] = copy_yuv(aux_image);
					free(aux_image);
				}
			}
			else{
				seq[f] = copy_yuv(video[1]);
			}

			printf(" --> Process completed\n");

			if(f == 0){
				f = bframes;

				back_ref_error = get_dec_err(dec, 1);
				video[0] = video[1];
			}
			else if(f == for_reference){
				if(f == back_reference + 1){
					f = frames;
				}
				else{
					f = back_reference + 1;

					for_ref_error = get_dec_err(dec, 1);

					video[2] = video[1];
				}
			}
			else if(back_reference < f && f < for_reference - 1){
				f++;

				safefree_yuv(&video[1]);
			}
			else if(f == for_reference - 1){
				if(f + 1 == frames - 1){
					f = frames;
				}
				else if(f + bframes + 1 < frames){
					f = f + bframes + 1;

					safefree((void **)&back_ref_error);
					back_ref_error = NULL;
					back_ref_error = for_ref_error;
					back_reference = for_reference;
					for_reference = f;
				}
				else{
					f = frames - 1;

					safefree((void **)&back_ref_error);
					back_ref_error = NULL;
					back_ref_error = for_ref_error;
					back_reference = for_reference;
					for_reference = f;
				}

				if(video[1] != video[0] && video[1] != video[2]){
					safefree_yuv(&video[1]);
				}
				safefree_yuv(&video[0]);
				video[0] = video[2];
			}

			free_decoder(dec);
		}


		if(video[1] != video[0] && video[1] != video[2]){
			if(video[1] != NULL) safefree_yuv(&video[1]);
		}
		if(video[0] != video[2]){
			safefree_yuv(&video[0]);
		}
		safefree_yuv(&video[2]);

		if(diff == 0){
			for(f = 0; f < frames; f++){
				write_yuv(seq[f], outfile);
				free(seq[f]->val);
				free(seq[f]);
			}
		}
		else{
			write_yuv(seq[0], outfile);
			img_diff = seq[0];
			for(f = 1; f < frames; f++){
				img_diff = sum_diff(img_diff, seq[f], f);
				write_yuv(img_diff, outfile);
				free(seq[f - 1]->val);
				free(seq[f - 1]);
			}
			free(seq[frames - 1]->val);
			free(seq[frames - 1]);

			safefree_yuv(&img_diff);
		}

		if(back_ref_error != for_ref_error){
			safefree((void **)&back_ref_error);
		}
		safefree((void **)&for_ref_error);

		free(seq);
	}
	else if(hevc == 1){
		IMAGE **seq = alloc_mem(frames * sizeof(IMAGE));
		int first_frame = 0;
		int **back_ref_error, **for_ref_error;
		int ***keep_error = (int ***)alloc_mem((bframes + 1) * sizeof(int **));

		for(f = 0; f < bframes + 1; f++) keep_error[f] = NULL;

		int (*bref)[5] = NULL;
		int conta = 0, final;

		final = ((frames - 1) % bframes != 0) ? bframes * ((int)((frames - 1) / bframes)) + 1 : frames;

		switch(bframes){
			case 3:
				bref = bref2;
				break;
			case 4:
				bref = bref3;
				break;
			case 5:
				bref = bref4;
				break;
			case 6:
				bref = bref5;
				break;
			case 7:
				bref = bref6;
				break;
			case 8:
				bref = bref7;
				break;
			case 9:
				bref = bref8;
				break;
			case 10:
				bref = bref9;
				break;
			default:
				bref = bref3;
				break;
		}

		f = 0;

		while(f < frames){
			printf("Decoding frame: %03d", f);

			if(f == 0){
				dec = init_decoder(fp, NULL, NULL, version, width, height, maxval, num_comp, num_group, prd_order[0], 0, 0, num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta);
			}
			else if((f - first_frame) % bframes == 0){
				dec = init_decoder(fp, back_ref_error, NULL, version, width, height, maxval, num_comp, num_group, prd_order[1], prd_order[2], 0, num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta);
			}
			else{
				dec = init_decoder(fp, back_ref_error, for_ref_error, version, width, height, maxval, num_comp, num_group, prd_order[3], prd_order[4], prd_order[5], num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta);
			}

			if (dec->f_huffman == 0){
				dec->rc = rc_init();
				rc_startdec(fp, dec->rc);
			}

			decode_class(fp, dec);
			decode_predictor(fp, dec);
			decode_threshold(fp, dec);

			dec->pmodels = init_pmodels(dec->num_group, dec->num_pmodel, dec->pm_accuracy, dec->pm_idx, dec->sigma, dec->maxval + 1);

			video[1] = decode_image(fp, video, dec);

			if(f > 0 && diff != 0){
				extra_info = decode_extra_info(fp, &num_pels);

				IMAGE *aux_image;
				if(num_pels != 0 && extra_info != NULL){
					aux_image = alloc_image(dec->width, dec->height, 255);
					int conta = 0;

					for(y = 0; y < dec->height; y++){
						for(x = 0; x < dec->width; x++){
							aux_image->val[y][x] = video[1]->val[y][x];

							if((video[1]->val[y][x] == 0 || video[1]->val[y][x] == 255) && conta < num_pels){
								aux_image->val[y][x] += extra_info[conta];
								conta++;
							}
						}
					}
				}

				if(extra_info != NULL) free(extra_info);

				if(num_pels == 0){
					seq[f] = copy_yuv(video[1]);
					safefree_yuv(&video[1]);
				}
				else{
					seq[f] = copy_yuv(aux_image);
					safefree_yuv(&aux_image);
				}
			}
			else{
				seq[f] = copy_yuv(video[1]);
				safefree_yuv(&video[1]);
			}

			printf(" --> Process completed\n");

			if(f == 0){
				f = bref[conta][0];

				keep_error[0] = get_dec_err(dec, 1);
				for_ref_error = NULL;
				video[2] = NULL;
				back_ref_error = keep_error[0];

				video[0] = seq[f + bref[conta][1]];

				conta++;
			}
			else{
				if(conta + 1 <= bframes){
					keep_error[f - first_frame] = get_dec_err(dec, 1);

					if(f + bref[conta][0] > 2 * first_frame && first_frame != 0){
						f = first_frame + bref[conta][0];
					}
					else{
						f = f + bref[conta][0];
					}

					video[0] = seq[f + bref[conta][1]];

					if(conta < bframes + 1) back_ref_error = keep_error[bref[conta][3]];

					if(bref[conta][2] != -1){
						video[2] = seq[f + bref[conta][2]];

						if(conta < bframes + 1) for_ref_error = keep_error[bref[conta][4]];
					}
					else{
						for_ref_error = NULL;
					}

					conta++;
				}
				else if(conta + 1 > bframes){
					if(first_frame + 2 * bframes >= frames){
						if(f + 1 < final && final != frames){
							int numbs = frames - final;
							int aux_bframes = bframes;
							bframes = numbs;

							conta = 1;

							video[0] = seq[f + 1];

							first_frame = f + 1;
							f = frames - 1;
							final = frames;

							switch(numbs){
								case 1:
									conta = bframes;

									video[0] = seq[frames - 2];

									first_frame = frames - 2;

									break;
								case 2:
									bref = bref1;
									break;
								case 3:
									bref = bref2;
									break;
								case 4:
									bref = bref3;
									break;
								case 5:
									bref = bref4;
									break;
								case 6:
									bref = bref5;
									break;
								case 7:
									bref = bref6;
									break;
								case 8:
									bref = bref7;
									break;
								case 9:
									bref = bref8;
									break;
								case 10:
									bref = bref9;
									break;
								default:
									bref = bref3;
									break;
							}

							video[2] = NULL;
							for_ref_error = NULL;

							back_ref_error = keep_error[aux_bframes];
							safefree((void **)&keep_error[0]);
							keep_error[0] = keep_error[aux_bframes];

							for(i = 1; i < aux_bframes; i++){
								safefree((void **)&keep_error[i]);
							}
						}
						else{
							break;
						}
					}
					else{
						conta = 1;
						f = first_frame + 2 * bframes;

						video[0] = seq[f + bref[0][1]];

						video[2] = NULL;
						for_ref_error = NULL;

						back_ref_error = keep_error[bframes];
						safefree((void **)&keep_error[0]);
						keep_error[0] = keep_error[bframes];

						for(i = 1; i < bframes; i++){
							safefree((void **)&keep_error[i]);
						}

						first_frame = first_frame + bframes;
					}
				}
			}

			free_decoder(dec);
			dec = NULL;
		}

		if(dec != NULL) free_decoder(dec);

		if(diff == 0){
			for(f = 0; f < frames; f++){
				write_yuv(seq[f], outfile);
				free(seq[f]->val);
				free(seq[f]);
			}
		}
		else{
			write_yuv(seq[0], outfile);
			img_diff = seq[0];
			for(f = 1; f < frames; f++){
				img_diff = sum_diff(img_diff, seq[f], f);
				write_yuv(img_diff, outfile);
				free(seq[f - 1]->val);
				free(seq[f - 1]);
			}
			free(seq[frames - 1]->val);
			free(seq[frames - 1]);

			safefree_yuv(&img_diff);
		}

		for(i = 0; i < bframes + 1; i++){
			safefree((void **)&keep_error[i]);
		}

		safefree((void **)&keep_error);

		free(seq);
	}

	fclose(fp);

	printf("\ncpu time: %.2f sec.\n", cpu_time());

	return(0);
}
