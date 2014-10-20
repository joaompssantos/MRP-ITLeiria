#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mrp.h"

extern POINT dyx[];
extern POINT idyx[];
extern double sigma_h[], sigma_a[];        
extern double qtree_prob[];

uint getbits(FILE *fp, int n){
	static int bitpos = 0;
	static uint bitbuf = 0;
	int x = 0;

	if (n <= 0) return (0);
	while (n > bitpos) {
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

DECODER *init_decoder(FILE *fp){
	DECODER *dec;
	int i, f;

	dec = (DECODER *)alloc_mem(sizeof(DECODER));

	if (getbits(fp, 16) != MAGIC_NUMBER) {
		fprintf(stderr, "Not a compressed file!\n");
		exit(1);
	}

	dec->version = getbits(fp, 8);
	dec->width = getbits(fp, 16);
	dec->height = getbits(fp, 16);
	dec->maxval = getbits(fp, 16);
	dec->frames = getbits(fp, 16);

	dec->num_class = (int *)malloc(dec->frames * sizeof(int));
	for(f = 0; f < dec->frames; f++){
		dec->num_class[f] = getbits(fp, 8);
	}

	dec->num_comp = getbits(fp, 4);
	dec->num_group = getbits(fp, 6);
	dec->prd_order = getbits(fp, 7);
	dec->inter_prd_order = getbits(fp, 8);
	dec->num_pmodel = getbits(fp, 6) + 1;
	dec->coef_precision = getbits(fp, 4) + 1;
	dec->pm_accuracy = getbits(fp, 3) - 1;
	dec->f_huffman = getbits(fp, 1);
	dec->quadtree_depth = (getbits(fp, 1))? QUADTREE_DEPTH : -1;

	dec->maxprd = dec->maxval << dec->coef_precision;

	dec->predictor = (int ***)malloc(dec->frames * sizeof(int **));
	i = dec->prd_order * 3 + 5;
	dec->predictor[0] = (int **) alloc_2d_array(dec->num_class[0], i, sizeof(int));
	i = (dec->prd_order + dec->inter_prd_order) * 3 + 5;
	for(f = 1; f < dec->frames; f++){
		dec->predictor[f] = (int **)alloc_2d_array(dec->num_class[f], i, sizeof(int));
	}

	dec->err = (int ***)alloc_3d_array(dec->height, dec->width, dec->frames, sizeof(int));

	dec->ctx_weight = (int **) alloc_mem(dec->frames * sizeof(int *));
	for(f = 0; f < dec->frames; f++){
		dec->ctx_weight[f] = init_ctx_weight();
	}

	// Quadtree map
	if (dec->quadtree_depth > 0) {
		int x, y, xx, yy;

		dec->qtmap = (char ****)alloc_2d_array(dec->frames, QUADTREE_DEPTH, sizeof(char**));

		for(f = 0; f < dec->frames; f++){
			yy = (dec->height + MAX_BSIZE - 1) / MAX_BSIZE;
			xx = (dec->width + MAX_BSIZE - 1) / MAX_BSIZE;

			for (i = dec->quadtree_depth - 1; i >= 0; i--) {
				dec->qtmap[f][i] = (char **)alloc_2d_array(yy, xx, sizeof(char));

				for (y = 0; y < yy; y++) {
					dec->qtmap[f][i][y] = (char *)malloc(xx * sizeof(char));
					for (x = 0; x < xx; x++) {
						dec->qtmap[f][i][y][x] = 0;
					}
				}

				yy <<= 1;
				xx <<= 1;
			}
		}
	}

	// Class and uquant arrays
	dec->class = (char ***)alloc_3d_array(dec->height, dec->width, dec->frames, sizeof(char));

	dec->uquant = (char ***)malloc(dec->frames * sizeof(char **));
	for(f = 0; f < dec->frames; f++){
		dec->uquant[f] = (char **)alloc_2d_array(dec->num_class[f], MAX_UPARA + 1, sizeof(char));
	}

	dec->pm_idx = (int **) alloc_mem(dec->frames * sizeof(int *));
	for(f = 0; f < dec->frames; f++){
		if (dec->num_pmodel > 1) {
			dec->pm_idx[f] = (int *)alloc_mem(dec->num_group * sizeof(int));
		}
		else {
			dec->pm_idx[f] = NULL;
		}
	}

	dec->spm = (PMODEL *) alloc_mem(dec->frames * sizeof(PMODEL));
	for(f = 0; f < dec->frames; f++){
		dec->spm[f].freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
		dec->spm[f].cumfreq = &(dec->spm[f].freq[MAX_SYMBOL]);
	}

	//Huffman coding
	if (dec->f_huffman == 1) {
		dec->sigma = sigma_h;
	}
	else {
		dec->sigma = sigma_a;
	}

	dec->mtfbuf = (int **)malloc(dec->frames * sizeof(int *));
	for(f = 0; f < dec->frames; f++){
		dec->mtfbuf[f] = (int *)malloc(dec->num_class[f] * sizeof(int));
	}

	dec->vlcs = (VLC ***) alloc_mem(dec->frames * sizeof(VLC **));
	dec->pmodels = (PMODEL ****) alloc_mem(dec->frames * sizeof(PMODEL ***));
	dec->rc = (RANGECODER **) alloc_mem(dec->frames * sizeof(RANGECODER *));

	return (dec);
}

int decode_vlc(FILE *fp, VLC *vlc){
	int i, k, min, off;
	uint code;

	code = min = off = k = 0;

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

int decode_golomb(FILE *fp, int m){
	int v = 0;

	while (getbits(fp, 1) == 0) {
		v++;
	}

	v = (v << m) | getbits(fp, m);

	return (v);
}

void decode_predictor(FILE *fp, DECODER *dec, int frame){
	int k, m, cl, coef, sgn, *coef_p, *scany_p, *scanx_p, *size_p;

	int prd_order = dec->prd_order;
	if(frame > 0) prd_order += dec->inter_prd_order;

	if (dec->f_huffman == 1) {
		for (k = 0; k < prd_order; k++) {
			m = getbits(fp, 4);

			for (cl = 0; cl < dec->num_class[frame]; cl++) {
				coef = decode_golomb(fp, m);

				if (coef > 0) {
					sgn = getbits(fp, 1);

					if (sgn) {
						coef = -coef;
					}
				}

				dec->predictor[frame][cl][k] = coef;
			}
		}
	}
	else {
		PMODEL *pm;

		pm = &dec->spm[frame];
		pm->size = MAX_COEF + 18;
		pm->cumfreq[MAX_COEF + 2] = 0;

		for(k = MAX_COEF + 2; k < pm->size; k++) {
			pm->freq[k] = 1;
			pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
		}

		for (k = 0; k < prd_order; k++) {
			m = rc_decode(fp, dec->rc[frame], pm, MAX_COEF + 2, MAX_COEF + 18) - (MAX_COEF + 2);
			set_spmodel(pm, MAX_COEF + 1, m);

			for (cl = 0; cl < dec->num_class[frame]; cl++) {
				coef = rc_decode(fp, dec->rc[frame], pm, 0, MAX_COEF + 1);

				if (coef > 0) {
					sgn = rc_decode(fp, dec->rc[frame], pm, MAX_COEF+2, MAX_COEF+4) - (MAX_COEF + 2);

					if (sgn) {
						coef = -coef;
					}
				}

				dec->predictor[frame][cl][k] = coef;
			}
		}
	}

	for (cl = 0; cl < dec->num_class[frame]; cl++) {
		coef_p = dec->predictor[frame][cl];
		scany_p = coef_p + prd_order;
		scanx_p = scany_p + prd_order;
		size_p = scanx_p + prd_order;
		size_p[1] = size_p[2] = size_p[3] = size_p[4] = 0;

		for (k = m = 0; k < dec->prd_order; k++) {
			if (coef_p[k] != 0) {
				coef_p[m] = coef_p[k];
				scany_p[m] = dyx[k].y;
				scanx_p[m] = dyx[k].x;

				if (dyx[k].y < size_p[1]) size_p[1] = dyx[k].y;
				if (dyx[k].x < size_p[2]) size_p[2] = dyx[k].x;
				if (dyx[k].x > size_p[3]) size_p[3] = dyx[k].x;

				m++;
			}
		}
		if(frame > 0){
			for(k = 0; k < dec->inter_prd_order; k++){
				coef_p[m] = coef_p[k + dec->prd_order];
				scany_p[m] = idyx[k].y;
				scanx_p[m] = idyx[k].x;

				if (idyx[k].y < size_p[1]) size_p[1] = idyx[k].y;
				if (idyx[k].x < size_p[2]) size_p[2] = idyx[k].x;
				if (idyx[k].x > size_p[3]) size_p[3] = idyx[k].x;
				if (idyx[k].y > size_p[4]) size_p[4] = idyx[k].y;

				m++;
			}

			m -= dec->inter_prd_order;
		}

		size_p[0] = m;
		size_p[1] = -size_p[1];
		size_p[2] = -size_p[2];
		size_p[3] = dec->width - size_p[3];
		size_p[4] = dec->height - size_p[4];
	}

	return;
}

void decode_threshold(FILE *fp, DECODER *dec, int frame){
	int cl, gr, m, u, k;

	if (dec->f_huffman == 1) {
		m = getbits(fp, 4);

		for (cl = 0; cl < dec->num_class[frame]; cl++) {
			k = u = 0;

			for (gr = 0; gr < dec->num_group; gr++) {
				if (k > MAX_UPARA || gr == dec->num_group - 1) {
					k = MAX_UPARA + 1;
				}
				else {
					if (getbits(fp, 1)) k += decode_golomb(fp, m) + 1;
				}

				for (; u < k; u++) dec->uquant[frame][cl][u] = gr;
			}
		}

		if (dec->num_pmodel > 1) {
			for (k = 1; (1 << k) < dec->num_pmodel; k++);

			for (gr = 0; gr < dec->num_group; gr++) {
				dec->pm_idx[frame][gr] = getbits(fp, k);
			}
		}
	}
	else {
		PMODEL *pm;

		pm = &dec->spm[frame];
		pm->size = 16;
		pm->cumfreq[0] = 0;

		for (k = 0; k < pm->size; k++) {
			pm->freq[k] = 1;
			pm->cumfreq[k + 1] = pm->cumfreq[k] + pm->freq[k];
		}

		m = rc_decode(fp, dec->rc[frame], pm, 0, pm->size);
		set_spmodel(pm, MAX_UPARA + 2, m);

		for (cl = 0; cl < dec->num_class[frame]; cl++) {
			k = u = 0;

			for (gr = 0; gr < dec->num_group; gr++) {
				if (k > MAX_UPARA || gr == dec->num_group - 1) {
					k = MAX_UPARA + 1;
				}
				else {
					k += rc_decode(fp, dec->rc[frame], pm, 0, pm->size - k);
				}

				for (; u < k; u++) dec->uquant[frame][cl][u] = gr;
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
				dec->pm_idx[frame][gr] = rc_decode(fp, dec->rc[frame], pm, 0, pm->size);
			}
		}
	}
	return;
}

void decode_qtindex(FILE *fp, DECODER *dec, int frame, VLC *vlc, PMODEL *cpm, int tly, int tlx, int blksize, int width, int level){
	int i, cl, y, x, bry, brx, ctx;
	char **qtmap;
	PMODEL *pm;

	brx = (tlx + blksize < dec->width) ? (tlx + blksize) : dec->width;
	bry = (tly + blksize < dec->height) ? (tly + blksize) : dec->height;

	if (tlx >= brx || tly >= bry) return;

	if (level > 0) {
		ctx = 0;
		qtmap = dec->qtmap[frame][level - 1];
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
			pm = &dec->spm[frame];
			i = rc_decode(fp, dec->rc[frame], pm, ctx, ctx + 2) - ctx;
		}

		if (i == 1) {
			qtmap[y][x] = 1;
			blksize >>= 1;
			decode_qtindex(fp, dec, frame, vlc, cpm, tly, tlx, blksize, width, level - 1);
			decode_qtindex(fp, dec, frame, vlc, cpm, tly, tlx + blksize, blksize, width, level - 1);
			decode_qtindex(fp, dec, frame, vlc, cpm, tly + blksize, tlx, blksize, width, level - 1);
			decode_qtindex(fp, dec, frame, vlc, cpm, tly + blksize, tlx + blksize, blksize, brx, level - 1);

			return;
		}
	}

	if (dec->f_huffman == 1) {
		i = decode_vlc(fp, vlc);
	}
	else {
		i = rc_decode(fp, dec->rc[frame], cpm, 0, cpm->size);
	}

	mtf_classlabel(dec->class[frame], dec->mtfbuf[frame], tly, tlx,blksize, width, dec->num_class[frame]);

	for (cl = 0; cl < dec->num_class[frame]; cl++) {
		if (dec->mtfbuf[frame][cl] == i) break;
	}

	for (y = tly; y < bry; y++) {
		for (x = tlx; x < brx; x++) {
			dec->class[frame][y][x] = cl;
		}
	}

	return;
}

void decode_class(FILE *fp, DECODER *dec, int frame){
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
		vlc->size = dec->num_class[frame];
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

		pm = &dec->spm[frame];
		if (dec->quadtree_depth > 0) {
			set_spmodel(pm, 7, -1);

			for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
				qtree_code[ctx] = rc_decode(fp, dec->rc[frame], pm, 0, pm->size);
			}
		}

		set_spmodel(pm, PMCLASS_LEVEL, -1);

		for (i = 0; i < dec->num_class[frame]; i++) {
			mtf_code[i] = rc_decode(fp, dec->rc[frame], pm, 0, pm->size);

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

		cpm->size = dec->num_class[frame];
		cpm->freq = (uint *)alloc_mem((cpm->size * 2 + 1) * sizeof(uint));
		cpm->cumfreq = &cpm->freq[cpm->size];
		cpm->cumfreq[0] = 0;

		for (i = 0; i < dec->num_class[frame]; i++) {
			p = exp(-log(2.0) * ((double)mtf_code[i] + 0.5) * PMCLASS_MAX/PMCLASS_LEVEL);
			cpm->freq[i] = p * (1 << 10);

			if (cpm->freq[i] <= 0) cpm->freq[i] = 1;

			cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
		}

		vlc = NULL;
	}

	for (i = 0; i < dec->num_class[frame]; i++) {
		dec->mtfbuf[frame][i] = i;
	}

	for (y = 0; y < dec->height; y += blksize) {
		for (x = 0; x < dec->width; x += blksize) {
			decode_qtindex(fp, dec, frame, vlc, cpm, y, x, blksize, dec->width, level);
		}
	}

	if (dec->f_huffman == 1) {
		free_vlc(vlc);
	}

	return;
}


int calc_udec(DECODER *dec, int frame, int y, int x){
	int rx, ry, k, u;
	int **err, *wt_p;

	u = 0;
	err = dec->err[frame];
	wt_p = dec->ctx_weight[frame];

	if (y >= UPEL_DIST && x >= UPEL_DIST && x <= dec->width - UPEL_DIST) {
		for (k = 0; k < NUM_UPELS; k++) {
			ry = y + dyx[k].y;
			rx = x + dyx[k].x;
			u += err[ry][rx] * (*wt_p++);
		}
	}
	else if (y == 0) {
		if (x == 0) {
			for (k = 0; k < NUM_UPELS; k++) {
				u += ((dec->maxval + 1) >> 2) * (*wt_p++);
			}
		}
		else {
			ry = 0;

			for (k =0; k < NUM_UPELS; k++) {
				rx = x + dyx[k].x;

				if (rx < 0) rx = 0;
				else if (rx >= x) rx = x - 1;

				u += err[ry][rx] * (*wt_p++);
			}
		}
	}
	else {
		if (x == 0) {
			for (k = 0; k < NUM_UPELS; k++) {
				ry = y + dyx[k].y;

				if (ry < 0) ry = 0;
				else if (ry >= y) ry = y - 1;

				rx = x + dyx[k].x;

				if (rx < 0) rx = 0;

				u += err[ry][rx] * (*wt_p++);
			}
		}
		else {
			for (k = 0; k < NUM_UPELS; k++) {
				ry = y + dyx[k].y;

				if (ry < 0) ry = 0;

				rx = x + dyx[k].x;

				if (rx < 0) rx = 0;
				else if (rx >= dec->width) rx = dec->width - 1;

				u += err[ry][rx] * (*wt_p++);
			}
		}
	}

	u >>= 6;

	if (u > MAX_UPARA) u = MAX_UPARA;

	return (u);
}

int calc_prd(IMAGE *img, DECODER *dec, int frame, int cl, int y, int x){
	int k, prd, prd_order, rx, ry, *coef_p, *scany_p, *scanx_p, *size_p;

	prd_order = dec->prd_order;
	if(frame > 0) prd_order += dec->inter_prd_order;

	coef_p = dec->predictor[frame][cl];
	scany_p = coef_p + prd_order;
	scanx_p = scany_p + prd_order;
	size_p = scanx_p + prd_order;
	prd_order = size_p[0];
	prd = 0;

	if (y >= size_p[1] && x >= size_p[2] && x < size_p[3]){
		for (k = 0; k < prd_order; k++) {
			ry = y + (*scany_p++);
			rx = x + (*scanx_p++);
			prd += (*coef_p++) * img->val[frame][ry][rx];
		}
	}
	else if (y == 0) {
		if (x == 0) {
			for (k = 0; k < prd_order; k++) {
				prd += *coef_p++;
			}

			prd *= ((img->maxval + 1) >> 1);
		}
		else {
			ry = 0;

			for (k = 0; k < prd_order; k++) {
				rx = x + (*scanx_p++);

				if (rx < 0) rx = 0;
				else if (rx >= x) rx = x - 1;

				prd += (*coef_p++) * img->val[frame][ry][rx];
			}
		}
	}
	else {
		if (x == 0) {
			for (k = 0; k < prd_order; k++) {
				ry = y + (*scany_p++);

				if (ry < 0) ry = 0;
				else if (ry >= y) ry = y - 1;

				rx = x + (*scanx_p++);

				if (rx < 0) rx = 0;

				prd += (*coef_p++) * img->val[frame][ry][rx];
			}
		}
		else {
			for (k = 0; k < prd_order; k++) {
				ry = y + (*scany_p++);

				if (ry < 0) ry = 0;

				rx = x + (*scanx_p++);

				if (rx < 0) rx = 0;
				else if (rx >= img->width) rx = img->width - 1;

				prd += (*coef_p++) * img->val[frame][ry][rx];
			}
		}
	}

	//Inter prediction calculation
	if(frame > 0 && dec->inter_prd_order > 0){
		if (y >= size_p[1] && x >= size_p[2] && x < size_p[3] && y < size_p[4]){
			prd += (*coef_p++) * img->val[frame - 1][y][x];

			for (k = 0; k < dec->inter_prd_order - 1; k++) {
				ry = y + (*scany_p++);
				rx = x + (*scanx_p++);

				prd += (*coef_p++) * img->val[frame - 1][ry][rx];
			}
		}
		else{
			prd += (*coef_p++) * img->val[frame - 1][y][x];

			for (k = 0; k < dec->inter_prd_order - 1; k++) {
				ry = y + (*scany_p++);
				rx = x + (*scanx_p++);

				if(ry < 0 || rx < 0 || ry >= dec->height || rx >= dec->width){
					prd += (*coef_p++) * img->val[frame - 1][y][x];
				}
				else{
					prd += (*coef_p++) * img->val[frame - 1][ry][rx];
				}
			}
		}
	}

	if (prd < 0) prd = 0;
	else if (prd > dec->maxprd) prd = dec->maxprd;

	return (prd);
}

IMAGE *decode_image(FILE *fp, IMAGE *img, DECODER *dec, int frame){
	int x, y, cl, gr, prd, u, e, E, p;

	if (dec->f_huffman == 1) {
		VLC *vlc;
		dec->vlcs[frame] = init_vlcs(dec->pmodels[frame], dec->num_group, 1);

		for (y = 0; y < dec->height; y++) {
			for (x = 0; x < dec->width; x++) {
				cl = dec->class[frame][y][x];
				u = calc_udec(dec, frame, y, x);
				gr = dec->uquant[frame][cl][u];
				prd = calc_prd(img, dec, frame, cl, y, x);
				prd >>= (dec->coef_precision - 1);
				p = (prd + 1) >> 1;
				vlc = &dec->vlcs[frame][gr][0];
				dec->err[frame][y][x] = E = decode_vlc(fp, vlc);
				e = E2e(E, p, prd & 1, dec->maxval);
				img->val[frame][y][x] = p + e;
			}
		}
	}
	else {
		PMODEL *pm;
		if (dec->pm_accuracy < 0) {
			for (y = 0; y < dec->height; y++) {
				for (x = 0; x < dec->width; x++) {
					cl = dec->class[frame][y][x];
					u = calc_udec(dec, frame, y, x);
					gr = dec->uquant[frame][cl][u];
					prd = calc_prd(img, dec, frame, cl, y, x);
					prd >>= (dec->coef_precision - 1);
					p = (prd + 1) >> 1;
					pm = dec->pmodels[frame][gr][0];
					dec->err[frame][y][x] = E = rc_decode(fp, dec->rc[frame], pm, 0, pm->size);
					e = E2e(E, p, prd & 1, dec->maxval);
					img->val[frame][y][x] = p + e;
				}
			}
		}
		else {
			int mask, shift, base;
			mask = (1 << dec->pm_accuracy) - 1;
			shift = dec->coef_precision - dec->pm_accuracy;

			for (y = 0; y < dec->height; y++) {
				for (x = 0; x < dec->width; x++) {
					cl = dec->class[frame][y][x];
					u = calc_udec(dec,frame, y, x);
					gr = dec->uquant[frame][cl][u];
					prd = calc_prd(img, dec, frame, cl, y, x);
					base = (dec->maxprd - prd + (1 << shift) / 2) >> shift;
					pm = dec->pmodels[frame][gr][0] + (base & mask);
					base >>= dec->pm_accuracy;
					p = rc_decode(fp, dec->rc[frame], pm, base, base+dec->maxval+1) - base;
					img->val[frame][y][x] = p;
					prd >>= (dec->coef_precision - 1);
					e = (p << 1) - prd;
					dec->err[frame][y][x] = (e > 0)? (e - 1) : (-e);
				}
			}
		}
	}
	return (img);
}

void write_yuv(IMAGE *img, char *filename){
	int i, j, f;
	FILE *fp;

	fp = fileopen(filename, "wb");

	for(f = 0; f < img->frames; f++){
		for (i = 0; i < img->height; i++) {
			for (j = 0; j < img->width; j++) {
				putc(img->val[f][i][j], fp);
			}
		}
		
		for (i = 0; i < img->height / 2; i++) {
			for (j = 0; j < img->width / 2; j++) {
				putc(128, fp);
			}
		}
		
		for (i = 0; i < img->height / 2; i++) {
			for (j = 0; j < img->width / 2; j++) {
				putc(128, fp);
			}
		}
	}

	fclose(fp);

	return;
}

int main(int argc, char **argv){
	int i, f;
	IMAGE *img;
	DECODER *dec;
	char *infile, *outfile;
	FILE *fp;

	cpu_time();
	setbuf(stdout, 0);
	infile = outfile = NULL;

	for (i = 1; i < argc; i++) {
		if (infile == NULL) {
			infile = argv[i];
		}
		else {
			outfile = argv[i];
		}
	}

	if (infile == NULL || outfile == NULL) {
		printf(BANNER"\n", 0.1 * VERSION);
		printf("usage: decmrp infile outfile\n");
		printf("infile:     Input file\n");
		printf("outfile:    Output file\n");
		exit(0);
	}

	fp = fileopen(infile, "rb");
	dec = init_decoder(fp);
	img = alloc_image(dec->width, dec->height, dec->frames, dec->maxval);

	printf("\nMRP-Video Decoder\n\n");
	// Print file characteristics to screen
	printf("%s -> %s (%dx%dx%d)\n", infile, outfile, img->width, img->height, img->frames);
	// Print coding parameters to screen
	printf("K = %d, J = %d, P = %d, V = %d, A = %d\n", dec->prd_order, dec->inter_prd_order, dec->coef_precision, dec->num_pmodel, dec->pm_accuracy);

	for(f = 0; f < dec->frames; f++){
		if (dec->f_huffman == 0) {
			dec->rc[f] = rc_init();
			rc_startdec(fp, dec->rc[f]);
		}

		decode_class(fp, dec, f);
		decode_predictor(fp, dec, f);
		decode_threshold(fp, dec, f);

		dec->pmodels[f] = init_pmodels(dec->num_group, dec->num_pmodel, dec->pm_accuracy, dec->pm_idx[f], dec->sigma, dec->maxval + 1);
		decode_image(fp, img, dec, f);
	}

	fclose(fp);

	write_yuv(img, outfile);

	printf("\nCPU time: %.2f sec.\n\n", cpu_time());
	return(0);
}
