#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include "mrp.h"

extern POINT dyx[];
extern POINT idyx[];
extern double sigma_h[], sigma_a[];        
extern double qtree_prob[];

int pixel_wise_predictor(ENCODER *enc){
	int f, i, j, aux;
	int resultado = 0;

	IMAGE *img = alloc_image(enc->width, enc->height, enc->frames, 255);

	for(f = 0; f < enc->frames; f++){
		for(i = 0; i < enc->height; i++){
			for(j = 0; j < enc->width; j++){
				if(f == 0){
					img->val[f][i][j] = enc->org[f][i][j];
				}
				else{
					aux = enc->org[f][i][j] + enc->org[f - 1][i][j];

					img->val[f][i][j] =  aux - 127;
				}
			}
		}
	}

	for(f = 1; f < enc->frames; f++){
		for(i = 0; i < enc->height; i++){
			for(j = 0; j < enc->width; j++){
				enc->org[f][i][j] = img->val[f][i][j];
			}
		}
	}

	free(img);
	return resultado;
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

DECODER *init_decoder(FILE *fp, int **error, int version, int width, int height, int maxval, int num_comp, int num_group, int prd_order, int inter_prd_order, int num_pmodel, int coef_precision, int pm_accuracy, int f_huffman, int quadtree_depth){
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
	dec->inter_prd_order = inter_prd_order;
	dec->num_pmodel = num_pmodel;
	dec->coef_precision = coef_precision;
	dec->pm_accuracy = pm_accuracy;
	dec->f_huffman = f_huffman;
	dec->quadtree_depth = quadtree_depth;
	dec->maxprd = dec->maxval << dec->coef_precision;

	dec->num_class = read_class(fp);

	dec->predictor = (int **)alloc_2d_array(dec->num_class, dec->prd_order + dec->inter_prd_order, sizeof(int));

	dec->err = (int ***)alloc_3d_array(dec->height, dec->width, 2, sizeof(int));

	if(error != NULL){
		int y, x;

		for(y = 0; y < dec->height; y++){
			for(x = 0; x < dec->width; x++){
				dec->err[0][y][x] = error[y][x];
			}
		}
	}

	dec->ctx_weight = init_ctx_weight(dec->prd_order, dec->inter_prd_order);

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

void read_header(FILE *fp, int *version, int *width, int *height, int *maxval, int *frames, int *num_comp, int *num_group, int *prd_order, int *intra_prd_order, int *inter_prd_order, int *num_pmodel, int *coef_precision, int *pm_accuracy, int *f_huffman, int *quadtree_depth){
	if (getbits(fp, 16) != MAGIC_NUMBER){
		fprintf(stderr, "Not a compressed file!\n");
		exit(1);
	}

	*version = getbits(fp, 8);
	*width = getbits(fp, 16);
	*height = getbits(fp, 16);
	*maxval = getbits(fp, 16);
	*frames = getbits(fp, 16);
	*num_comp = getbits(fp, 4);
	*num_group = getbits(fp, 6);
	*prd_order = getbits(fp, 7);
	*intra_prd_order = getbits(fp, 8);
	*inter_prd_order = getbits(fp, 8);
	*num_pmodel = getbits(fp, 6) + 1;
	*coef_precision = getbits(fp, 4) + 1;
	*pm_accuracy = getbits(fp, 3) - 1;
	*f_huffman = getbits(fp, 1);
	*quadtree_depth = (getbits(fp, 1))? QUADTREE_DEPTH : -1;
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

void decode_predictor(FILE *fp, DECODER *dec, int f){
	int k, m, cl, coef, sgn;

	int prd_order = dec->prd_order + dec->inter_prd_order;

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
	int rx, ry, k, u = 0;

	int order = NUM_UPELS;

	if(dec->prd_order < NUM_UPELS && dec->inter_prd_order != 0){
		order = dec->prd_order;
	}

	if (y >= UPEL_DIST && x >= UPEL_DIST && x <= dec->width - UPEL_DIST){
		for (k = 0; k < order; k++){
			ry = y + dyx[k].y;
			rx = x + dyx[k].x;
			u += dec->err[1][ry][rx] * dec->ctx_weight[k];
		}
	}
	else if (y == 0){
		if (x == 0){
			for (k = 0; k < NUM_UPELS; k++){
				u += ((dec->maxval + 1) >> 2) * dec->ctx_weight[k];
			}
		}
		else{
			ry = 0;

			for (k =0; k < order; k++){
				rx = x + dyx[k].x;

				if (rx < 0) rx = 0;
				else if (rx >= x) rx = x - 1;

				u += dec->err[1][ry][rx] * dec->ctx_weight[k];
			}
		}
	}
	else{
		if (x == 0){
			for (k = 0; k < order; k++){
				ry = y + dyx[k].y;

				if (ry < 0) ry = 0;
				else if (ry >= y) ry = y - 1;

				rx = x + dyx[k].x;

				if (rx < 0) rx = 0;

				u += dec->err[1][ry][rx] * dec->ctx_weight[k];
			}
		}
		else{
			for (k = 0; k < order; k++){
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
	if(dec->prd_order < NUM_UPELS && dec->inter_prd_order != 0){
		int iorder = NUM_UPELS - dec->prd_order;

		u += dec->err[0][y][x] * dec->ctx_weight[order];

		if (y >= UPEL_DIST && x >= UPEL_DIST && x < dec->width - UPEL_DIST && y < dec->height - UPEL_DIST){
			for (k = 0; k < iorder - 1; k++){
				ry = y + idyx[k].y;
				rx = x + idyx[k].x;

				u += dec->err[0][ry][rx] * dec->ctx_weight[order + k + 1];
			}
		}
		else{
			for (k = 0; k < iorder - 1; k++){
				ry = y + idyx[k].y;
				rx = x + idyx[k].x;

				if(ry < 0 || rx < 0 || ry >= dec->height || rx >= dec->width){
					u += dec->err[0][y][x] * dec->ctx_weight[order + k + 1];
				}
				else{
					u += dec->err[0][ry][rx] * dec->ctx_weight[order + k + 1];
				}
			}
		}
	}

	u >>= 6;

	if (u > MAX_UPARA) u = MAX_UPARA;

	return (u);
}

int calc_prd(IMAGE *video[2], DECODER *dec, int cl, int y, int x){
	int k, prd, rx, ry;

	int dy, dx, min_dx, max_dx, min_dy;
	int imin_dx, imax_dx, imin_dy, imax_dy;
	int min_abs_dx, max_abs_dx, min_abs_dy, max_abs_dy;

	min_dx = max_dx = min_dy = 0;
	imin_dx = imax_dx = imin_dy = imax_dy = 0;

	//Values to check for special cases
	for (k = 0; k < dec->prd_order; k++){
		dy = dyx[k].y;
		dx = dyx[k].x;

		if (dy < min_dy) min_dy = dy;
		if (dx < min_dx) min_dx = dx;
		if (dx > max_dx) max_dx = dx;
	}

	for (k = 0; k < dec->inter_prd_order - 1; k++){
		dy = idyx[k].y;
		dx = idyx[k].x;

		if (dy < imin_dy) imin_dy = dy;
		if (dy > imax_dy) imax_dy = dy;
		if (dx < imin_dx) imin_dx = dx;
		if (dx > imax_dx) imax_dx = dx;
	}

	max_abs_dy = dec->height - imax_dy;
	if(dec->inter_prd_order <= 1) max_abs_dy = dec->height;

	if(min_dy < imin_dy){
		min_abs_dy = -min_dy;
	}
	else{
		min_abs_dy = -imin_dy;
	}

	if(max_dx < imax_dx){
		max_abs_dx = dec->width - imax_dx;
	}
	else{
		max_abs_dx = dec->width - max_dx;
	}
	if(min_dx < imin_dx){
		min_abs_dx = -min_dx;
	}
	else{
		min_abs_dx = -imin_dx;
	}

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

	//Inter prediction calculation
	if(dec->inter_prd_order == 1){
		prd += dec->predictor[cl][dec->prd_order] * video[0]->val[y][x];
	}
	if(dec->inter_prd_order > 1){
		prd += dec->predictor[cl][dec->prd_order] * video[0]->val[y][x];

		if (y >= min_abs_dy && x >= min_abs_dx && x < max_abs_dx && y < max_abs_dy){
			for (k = 0; k < dec->inter_prd_order - 1; k++){
				ry = y + idyx[k].y;
				rx = x + idyx[k].x;

				prd += dec->predictor[cl][k + dec->prd_order + 1] * video[0]->val[ry][rx];
			}
		}
		else{
			for (k = 0; k < dec->inter_prd_order - 1; k++){
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

	if (prd < 0) prd = 0;
	else if (prd > dec->maxprd) prd = dec->maxprd;

	return (prd);
}

IMAGE *decode_image(FILE *fp, IMAGE *video[2], DECODER *dec){
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

void write_yuv(IMAGE *img, char *filename){
	int i, j;
	FILE *fp;

	fp = fileopen(filename, "ab");

	for (i = 0; i < img->height; i++){
		for (j = 0; j < img->width; j++){
			putc(img->val[i][j], fp);
		}
	}

	for (i = 0; i < img->height / 2; i++){
		for (j = 0; j < img->width / 2; j++){
			putc(128, fp);
		}
	}

	for (i = 0; i < img->height / 2; i++){
		for (j = 0; j < img->width / 2; j++){
			putc(128, fp);
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

int main(int argc, char **argv){
	int i, f, **error;
	int version, width, height, maxval, frames, num_comp, num_group, prd_order, intra_prd_order, inter_prd_order, num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth;
	IMAGE *video[2];
	DECODER *dec;
	char *infile, *outfile;
	FILE *fp;

	cpu_time();
	setbuf(stdout, 0);
	infile = outfile = NULL;

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
	read_header(fp, &version, &width, &height, &maxval, &frames, &num_comp, &num_group, &prd_order, &intra_prd_order, &inter_prd_order, &num_pmodel, &coef_precision, &pm_accuracy, &f_huffman, &quadtree_depth);

	printf("\nMRP-Video Decoder\n\n");
	// Print file characteristics to screen
	printf("%s -> %s (%dx%dx%d)\n", infile, outfile, width, height, frames);
	// Print coding parameters to screen
	printf("K = %d, L = %d, J = %d, P = %d, V = %d, A = %d\n", prd_order, intra_prd_order, inter_prd_order, coef_precision, num_pmodel, pm_accuracy);

	for(f = 0; f < frames; f++){
		printf("Decoding frame: %03d", f);

		if(f == 0){
			dec = init_decoder(fp, NULL, version, width, height, maxval, num_comp, num_group, prd_order, 0, num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth);
		}
		else{
			video[0] = read_yuv(outfile, height, width, f - 1);

			dec = init_decoder(fp, error, version, width, height, maxval, num_comp, num_group, intra_prd_order, inter_prd_order, num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth);
			free(error);
		}

		if (dec->f_huffman == 0){
			dec->rc = rc_init();
			rc_startdec(fp, dec->rc);
		}

		decode_class(fp, dec);
		decode_predictor(fp, dec, f);
		decode_threshold(fp, dec);

		dec->pmodels = init_pmodels(dec->num_group, dec->num_pmodel, dec->pm_accuracy, dec->pm_idx, dec->sigma, dec->maxval + 1);

		video[1] = decode_image(fp, video, dec);

		write_yuv(video[1], outfile);

		free(video[1]->val);
		free(video[1]);
		if(f > 0){
			free(video[0]->val);
			free(video[0]);
		}

		error = get_dec_err(dec, 1);
		free_decoder(dec);

		printf(" --> Process completed\n");
	}

	free(error);
	fclose(fp);

	printf("\ncpu time: %.2f sec.\n", cpu_time());

	return(0);
}
