/*
“Commons Clause” License Condition v1.0

The Software is provided to you by the Licensor under the License, as defined below, subject to the following condition.

Without limiting other conditions in the License, the grant of rights under the License will not include, and the License does not grant to you, the right to Sell the Software.

For purposes of the foregoing, “Sell” means practicing any or all of the rights granted to you under the License to provide to third parties, for a fee or other consideration (including without limitation fees for hosting or consulting/ support services related to the Software), a product or service whose value derives, entirely or substantially, from the functionality of the Software. Any license notice or attribution required by the License must also include this Commons Clause License Condition notice.

Software: Hierarchical Minimum Rate Predictors

License: BSD-3-Clause with Commons Clause

Licensor: Instituto de Telecomunicações

Copyright © 2023 Instituto de Telecomunicações. All Rights Reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include <sys/stat.h>
#include <unistd.h>
#include "mrp.h"
#include "mrp_config.h"

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
    return((int) getbits(fp, 8));
}

void read_hilevel(int *hilevel, int *frame, FILE *fp) {
    *hilevel = (int) getbits(fp, 8);
    *frame = (int) getbits(fp, 8);
}

void decode_disparity(FILE *fp, DECODER *dec, int disp_blk_size, int max_disparity) {
    int count = 0;
    int ref_vectors[2][3];
    int num_values = 4 * max_disparity + 1;

    PMODEL *pm = &dec->spm;
    // Uniform initialization of the probability model
    set_spmodel(pm, num_values, -1);

    for (int r = 0; r < dec->no_refsai; r++) {
        dec->disparity_vectors[r] = (int ***) alloc_3d_array((int) ceil(dec->ts[HEIGHT] * 1.0 / disp_blk_size), (int) ceil(dec->ts[WIDTH] * 1.0 / disp_blk_size), 2, sizeof(int));

        for (int t = 0; t < ceil(dec->ts[HEIGHT] * 1.0 / disp_blk_size); t++) {
            for (int s = 0; s < ceil(dec->ts[WIDTH] * 1.0 / disp_blk_size); s++) {
                count = 0;

                if (s > 0) {
                    ref_vectors[0][count] = dec->disparity_vectors[r][t][s - 1][0];
                    ref_vectors[1][count] = dec->disparity_vectors[r][t][s - 1][1];

                    count++;
                }

                if (t > 0) {
                    ref_vectors[0][count] = dec->disparity_vectors[r][t - 1][s][0];
                    ref_vectors[1][count] = dec->disparity_vectors[r][t - 1][s][1];

                    count++;
                }

                if (t > 1 && s < ceil(dec->ts[WIDTH] * 1.0 / disp_blk_size)) {
                    ref_vectors[0][count] = dec->disparity_vectors[r][t - 1][s + 1][0];
                    ref_vectors[1][count] = dec->disparity_vectors[r][t - 1][s + 1][1];

                    count++;
                }

                dec->disparity_vectors[r][t][s][0] = rc_decode(fp, dec->rc, pm, 0, num_values);
                update_pmodel(pm, dec->disparity_vectors[r][t][s][0], num_values);

                dec->disparity_vectors[r][t][s][1] = rc_decode(fp, dec->rc, pm, 0, num_values);
                update_pmodel(pm, dec->disparity_vectors[r][t][s][1], num_values);

                if (count > 0) {
                    dec->disparity_vectors[r][t][s][0] = dec->disparity_vectors[r][t][s][0] + floor(median(ref_vectors[0], count)) - 2 * max_disparity;
                    dec->disparity_vectors[r][t][s][1] = dec->disparity_vectors[r][t][s][1] + floor(median(ref_vectors[1], count)) - 2 * max_disparity;
                }
                else {
                    dec->disparity_vectors[r][t][s][0] = dec->disparity_vectors[r][t][s][0] - 2 * max_disparity;
                    dec->disparity_vectors[r][t][s][1] = dec->disparity_vectors[r][t][s][1] - 2 * max_disparity;
                }
            }
        }
    }
}

DECODER *init_decoder(FILE *fp, const int *ts, int maxval, int num_class, int num_comp, int num_group, int prd_order,
                      int sai_prd_order, int coef_precision, int f_huffman, int quadtree_depth, int num_pmodel,
                      int pm_accuracy, int delta, int depth, int frame, LF4D *lf, LF4D *save_err, int no_refsai,
                      int *reference_list, char *debug_path, bool use_disp, int disp_blk_size, int max_disparity) {
    DECODER *dec;
    int i;
    int sai_coord[2] = {0, 0};
    frame2coordinates(sai_coord, frame, lf->u);

    dec = (DECODER *) alloc_mem(sizeof(DECODER));

    dec->ts[HEIGHT] = ts[HEIGHT];
    dec->ts[WIDTH] = ts[WIDTH];

    dec->maxval = maxval;
    dec->num_comp = num_comp;
    dec->num_group = num_group;
    dec->no_refsai = no_refsai;

    dec->prd_order = prd_order;
    dec->full_prd_order = dec->prd_order;
    for (i = 0; i < dec->no_refsai; i++) {
        dec->sai_prd_order = sai_prd_order;
        dec->full_prd_order += dec->sai_prd_order;
    }

    dec->num_pmodel = num_pmodel;
    dec->coef_precision = coef_precision;
    dec->pm_accuracy = pm_accuracy;
    dec->f_huffman = f_huffman;
    dec->quadtree_depth = quadtree_depth;
    dec->maxprd = dec->maxval << dec->coef_precision;
    dec->delta = delta;
    dec->depth = depth;

    if (debug_path != NULL) {
        dec->debug_path = (char *) alloc_mem((strlen(debug_path) + 1) * sizeof(char));
        strcpy(dec->debug_path, debug_path);
    }
    else {
        dec->debug_path = NULL;
    }

    dec->num_class = num_class;

    dec->predictor = (int **) alloc_2d_array(dec->num_class, dec->full_prd_order, sizeof(int));

    // Memory allocation for the original image
    dec->org = (img_t ***) alloc_3d_array(dec->no_refsai + 1, dec->ts[HEIGHT] + 1, dec->ts[WIDTH], sizeof(img_t));
    // Memory allocation for the error structure
    dec->err = (int ***) alloc_3d_array(dec->no_refsai + 1, dec->ts[HEIGHT] + 1, dec->ts[WIDTH], sizeof(int));

    // Auxiliary values
    dec->org[0][dec->ts[HEIGHT]][0] = (dec->maxval + 1) >> 1;
    dec->err[0][dec->ts[HEIGHT]][0] = (dec->maxval + 1) >> 2;

    // Copy reference images to encoder struct, initialize group array
    for (i = 0; i < dec->no_refsai; i++) {
        frame2coordinates(sai_coord, reference_list[i], lf->u);

        for (int t = 0; t < dec->ts[HEIGHT]; t++) {
            for (int s = 0; s < dec->ts[WIDTH]; s++) {

                dec->org[i + 1][t][s] = lf->val[sai_coord[0]][sai_coord[1]][t][s];
                dec->err[i + 1][t][s] = save_err->val[sai_coord[0]][sai_coord[1]][t][s];
            }
        }

        // Auxiliary values
        dec->org[i][dec->ts[HEIGHT]][0] = (dec->maxval + 1) >> 1;
        dec->err[i][dec->ts[HEIGHT]][0] = (dec->maxval + 1) >> 2;

        // Initialization
        dec->disparity_vectors[i] = NULL;
    }

    //Initiation of the reference offset
    dec->intra_roff = init_ref_offset(dec->ts, 0, dec->prd_order, NULL, 0);

    // Keeps the weights used for the residue encoding context
    dec->intra_ctx_weight = init_ctx_weight(0, dec->prd_order, dec->delta);

    if (dec->f_huffman == 0) {
        dec->rc = rc_init();
        rc_startdec(fp, dec->rc);
    }

    dec->spm.freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
    dec->spm.cumfreq = &(dec->spm.freq[MAX_SYMBOL]);

    if (dec->no_refsai > 0 && use_disp) {
        decode_disparity(fp, dec, disp_blk_size, max_disparity);
    }

    for (i = 0; i < dec->no_refsai; i++) {
        if (use_disp) {
            dec->sai_ref_roff[i] = init_ref_offset(dec->ts, i + 1, dec->sai_prd_order, dec->disparity_vectors[i], disp_blk_size);
        }
        else{
            dec->sai_ref_roff[i] = init_ref_offset(dec->ts, i + 1, dec->sai_prd_order, NULL, 0);
        }

        dec->sai_ref_ctx_weight[i] = init_ctx_weight(i + 1, dec->sai_prd_order, dec->delta);
    }

    // Quadtree map
    if (dec->quadtree_depth > 0) {
        int t, s, tt, ss;

        tt = (dec->ts[HEIGHT] + MAX_BSIZE - 1) / MAX_BSIZE;
        ss = (dec->ts[WIDTH] + MAX_BSIZE - 1) / MAX_BSIZE;

        for (i = dec->quadtree_depth - 1; i >= 0; i--) {
            dec->qtmap[i] = (char **) alloc_2d_array(tt, ss, sizeof(char));

            for (t = 0; t < tt; t++) {
                for (s = 0; s < ss; s++) {
                    dec->qtmap[i][t][s] = 0;
                }
            }

            tt <<= 1;
            ss <<= 1;
        }
    }

    // Class and uquant arrays
    dec->class = (char **) alloc_2d_array(dec->ts[HEIGHT], dec->ts[WIDTH], sizeof(char));
    dec->uquant = (char **) alloc_2d_array(dec->num_class, MAX_UPARA + 1, sizeof(char));

    if (dec->num_pmodel > 1) {
        dec->pm_idx = (int *) alloc_mem(dec->num_group * sizeof(int));
    }
    else {
        dec->pm_idx = NULL;
    }

    //Huffman coding
    if (dec->f_huffman == 1) {
        dec->sigma = sigma_h;
    }
    else {
        dec->sigma = sigma_a;
    }

    dec->mtfbuf = (int *) alloc_mem(dec->num_class * sizeof(int));

    return (dec);
}

void free_decoder(DECODER *dec) {
    int i, j, gr, num_subpm;

    free(dec->predictor);
    free(dec->org);
    free(dec->err);

    for (i = 0; i < dec->no_refsai; i++) {
        if (dec->disparity_vectors[i] != NULL) free(dec->disparity_vectors[i]);
    }

    free(dec->intra_ctx_weight);
    free_ref_offset(dec->ts, 0, dec->prd_order, dec->intra_roff);

    for (i = 0; i < dec->no_refsai; i++) {
        free(dec->sai_ref_ctx_weight[i]);
        free_ref_offset(dec->ts, i + 1, dec->sai_prd_order, dec->sai_ref_roff[i]);
    }

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

    free(dec->debug_path);
    free(dec);
}

void read_header(FILE *fp, int vu[2], int ts[2], int *maxval, int *depth, int *num_comp, int *num_group, int *prd_order,
                 int *sai_prd_order, int *num_pmodel, int *coef_precision, int *pm_accuracy, int *f_huffman,
                 int *quadtree_depth, int *delta, int *hist_bytes, int *hilevels, int *no_ra_regions, bool *use_current,
                 bool *use_disp, int *disp_blk_size, int *max_disparity, int *sai_references, int *sai_distance_threshold) {
    if (getbits(fp, 16) != MAGIC_NUMBER) {
        fprintf(stderr, "Not a compressed file!\n");
        exit(1);
    }

    vu[HEIGHT] = (int) getbits(fp, 16);
    vu[WIDTH]  = (int) getbits(fp, 16);
    ts[HEIGHT] = (int) getbits(fp, 16);
    ts[WIDTH]  = (int) getbits(fp, 16);
    *maxval    = (int) getbits(fp, 16);
    *depth     = (int) getbits(fp, 6);
    *num_comp  = (int) getbits(fp, 2);
    *num_group = (int) getbits(fp, 6);
    *hilevels  = (int) getbits(fp, 8);

    *no_ra_regions = (int) getbits(fp, 8);
    *use_current = (bool) getbits(fp, 1);
    *use_disp = (bool) getbits(fp, 1);
    if (*use_disp == 1) {
        *disp_blk_size = (int) getbits(fp, 8);
        *max_disparity = (int) getbits(fp, 8);
    }
    *prd_order      = (int) getbits(fp, 8);
    *sai_prd_order  = (int) getbits(fp, 8);

    *sai_references         = (int) getbits(fp, 4);
    *sai_distance_threshold = (int) getbits(fp, 4);

    *delta          = (int) getbits(fp, 8);
    *num_pmodel     = (int) getbits(fp, 6) + 1;
    *coef_precision = (int) getbits(fp, 4) + 1;
    *pm_accuracy    = (int) getbits(fp, 4) - 1;
    *f_huffman      = (int) getbits(fp, 1);
    *quadtree_depth = (getbits(fp, 1)) ? QUADTREE_DEPTH : -1;
    *hist_bytes     = (int) getbits(fp, 16);
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
                    sgn = rc_decode(fp, dec->rc, pm, MAX_COEF + 2, MAX_COEF + 4) - (MAX_COEF + 2);

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

void decode_qtindex(FILE *fp, DECODER *dec, VLC *vlc, PMODEL *cpm, int tlt, int tls, int blksize, int width, int level) {
    int i, cl, t, s, brt, brs, ctx;
    char **qtmap;
    PMODEL *pm;

    brt = (tlt + blksize < dec->ts[HEIGHT]) ? (tlt + blksize) : dec->ts[HEIGHT];
    brs = (tls + blksize < dec->ts[WIDTH]) ? (tls + blksize) : dec->ts[WIDTH];

    if (tlt >= brt || tls >= brs) return;

    if (level > 0) {
        ctx = 0;
        qtmap = dec->qtmap[level - 1];

        t = (tlt / MIN_BSIZE) >> level;
        s = (tls / MIN_BSIZE) >> level;

        if (t > 0) {
            if (qtmap[t - 1][s] == 1) ctx++;
            if (brs < width && qtmap[t - 1][s + 1] == 1) ctx++;
        }

        if (s > 0 && qtmap[t][s - 1] == 1) ctx++;

        ctx = ((level - 1) * 4 + ctx) << 1;

        if (dec->f_huffman == 1) {
            i = getbits(fp, 1);
        }
        else {
            pm = &dec->spm;
            i = rc_decode(fp, dec->rc, pm, ctx, ctx + 2) - ctx;
        }

        if (i == 1) {
            qtmap[t][s] = 1;

            blksize >>= 1;

            // v and u
            decode_qtindex(fp, dec, vlc, cpm, tlt, tls, blksize, width, level - 1);
            decode_qtindex(fp, dec, vlc, cpm, tlt, tls + blksize, blksize, width, level - 1);
            decode_qtindex(fp, dec, vlc, cpm, tlt + blksize, tls, blksize, width, level - 1);
            decode_qtindex(fp, dec, vlc, cpm, tlt + blksize, tls + blksize, blksize, brs, level - 1);

            return;
        }
    }

    if (dec->f_huffman == 1) {
        i = decode_vlc(fp, vlc);
    }
    else {
        i = rc_decode(fp, dec->rc, cpm, 0, cpm->size);
    }

    mtf_classlabel(dec->class, dec->mtfbuf, tlt, tls, blksize, width, dec->num_class);

    for (cl = 0; cl < dec->num_class; cl++) {
        if (dec->mtfbuf[cl] == i) break;
    }

    for (t = tlt; t < brt; t++) {
        for (s = tls; s < brs; s++) {
            dec->class[t][s] = (char) cl;
        }
    }
}

void decode_class(FILE *fp, DECODER *dec) {
    int i, j, t, s, blksize, level;
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
        vlc = (VLC *) alloc_mem(sizeof(VLC));
        vlc->size = dec->num_class;
        vlc->max_len = 16;
        vlc->len = (int *) alloc_mem(vlc->size * sizeof(int));

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
        cpm->freq = (uint *) alloc_mem((cpm->size * 2 + 1) * sizeof(uint));
        cpm->cumfreq = &cpm->freq[cpm->size];
        cpm->cumfreq[0] = 0;

        for (i = 0; i < dec->num_class; i++) {
            p = exp(-log(2.0) * ((double) mtf_code[i] + 0.5) * PMCLASS_MAX / PMCLASS_LEVEL);
            cpm->freq[i] = (uint) (p * (1 << 10));

            if (cpm->freq[i] <= 0) cpm->freq[i] = 1;

            cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
        }

        vlc = NULL;
    }

    for (i = 0; i < dec->num_class; i++) {
        dec->mtfbuf[i] = i;
    }

    for (t = 0; t < dec->ts[HEIGHT]; t += blksize) {
        for (s = 0; s < dec->ts[WIDTH]; s += blksize) {
            decode_qtindex(fp, dec, vlc, cpm, t, s, blksize, dec->ts[WIDTH], level);
        }
    }

    if (dec->f_huffman == 1) {
        free_vlc(vlc);
    }
    else {
        free(cpm->freq);
    }
}

int calc_udec(DECODER *dec, int t, int s) {
    int uu, k, r, *err_p;
    int *intra_wt_p;
    int *intra_roff_p = NULL;
    int *sai_wt_p[dec->no_refsai];
    int *sai_roff_p[dec->no_refsai];

    err_p = &dec->err[0][t][s];

    intra_roff_p = dec->intra_roff[t][s];
    intra_wt_p = dec->intra_ctx_weight;

    for (r = 0; r < dec->no_refsai; r++) {
        sai_roff_p[r] = dec->sai_ref_roff[r][t][s];
        sai_wt_p[r] = dec->sai_ref_ctx_weight[r];
    }

    uu = 0;

    for (k = 0; k < dec->prd_order; k++) {
        uu += err_p[*intra_roff_p++] * (*intra_wt_p++);
    }

    for (r = 0; r< dec->no_refsai; r++) {
        for (k = 0; k < dec->sai_prd_order; k++) {
            uu += err_p[*sai_roff_p[r]++] * (*sai_wt_p[r]++);
        }
    }

    uu >>= 6;

    if (uu > MAX_UPARA) uu = MAX_UPARA;

    return (uu);
}

int calc_prd(DECODER *dec, int cl, int t, int s) {
    int k, r, prd, *coef_p;
    img_t *org_p;
    int *intra_roff_p = NULL;
    int *sai_roff_p[dec->no_refsai];

    intra_roff_p = dec->intra_roff[t][s];

    for (r = 0; r < dec->no_refsai; r++) {
        sai_roff_p[r] = dec->sai_ref_roff[r][t][s];
    }

    org_p = &dec->org[0][t][s];
    coef_p = dec->predictor[cl];

    prd = 0;

    for (k = 0; k < dec->prd_order; k++) {
        prd += org_p[*intra_roff_p++] * (*coef_p++);
    }

    for (r = 0; r < dec->no_refsai; r++) {
        for (k = 0; k < dec->sai_prd_order; k++) {
            prd += org_p[*sai_roff_p[r]++] * (*coef_p++);
        }
    }

    if (prd < 0) prd = 0;
    else if (prd > dec->maxprd) prd = dec->maxprd;

    return (prd);
}

void decode_image(FILE *fp, DECODER *dec) {
    int t, s, cl, gr, prd, uu, e, E, p;

    if (dec->f_huffman == 1) {
        VLC *vlc;
        dec->vlcs = init_vlcs(dec->pmodels, dec->num_group, 1);

        for (t = 0; t < dec->ts[HEIGHT]; t++) {
            for (s = 0; s < dec->ts[WIDTH]; s++) {
                cl = dec->class[t][s];
                uu = calc_udec(dec, t, s);
                gr = dec->uquant[cl][uu];
                prd = calc_prd(dec, cl, t, s);
                prd >>= (dec->coef_precision - 1);
                p = (prd + 1) >> 1;
                vlc = &dec->vlcs[gr][0];
                dec->err[0][t][s] = E = decode_vlc(fp, vlc);
                e = E2e(E, p, prd & 1, dec->maxval);
                dec->org[0][t][s] = (img_t) (p + e);
            }
        }
    }
    else {
        PMODEL *pm;
        if (dec->pm_accuracy < 0) {
            for (t = 0; t < dec->ts[HEIGHT]; t++) {
                for (s = 0; s < dec->ts[WIDTH]; s++) {
                    cl = dec->class[t][s];
                    uu = calc_udec(dec, t, s);
                    gr = dec->uquant[cl][uu];
                    prd = calc_prd(dec, cl, t, s);
                    prd >>= (dec->coef_precision - 1);
                    p = (prd + 1) >> 1;
                    pm = dec->pmodels[gr][0];
                    dec->err[0][t][s] = E = rc_decode(fp, dec->rc, pm, 0, pm->size);
                    e = E2e(E, p, prd & 1, dec->maxval);
                    dec->org[0][t][s] = (img_t) (p + e);
                }
            }
        }
        else {
            int mask, shift, base;
            mask = (1 << dec->pm_accuracy) - 1;
            shift = dec->coef_precision - dec->pm_accuracy;

            for (t = 0; t < dec->ts[HEIGHT]; t++) {
                for (s = 0; s < dec->ts[WIDTH]; s++) {
                    cl = dec->class[t][s];
                    uu = calc_udec(dec, t, s);
                    gr = dec->uquant[cl][uu];
                    prd = calc_prd(dec, cl, t, s);
                    base = (dec->maxprd - prd + (1 << shift) / 2) >> shift;
                    pm = dec->pmodels[gr][0] + (base & mask);
                    base >>= dec->pm_accuracy;
                    p = rc_decode(fp, dec->rc, pm, base, base + dec->maxval + 1) - base;
                    dec->org[0][t][s] = (img_t) p;
                    prd >>= (dec->coef_precision - 1);
                    e = (p << 1) - prd;
                    dec->err[0][t][s] = (e > 0) ? (e - 1) : (-e);
                }
            }
        }
    }
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
LF4D *histogram_unpacking(LF4D *lf, const int *backward_table) {
    int v, u, t, s;
    LF4D *u_lf = alloc_lf4d(lf->v, lf->u, lf->t, lf->s, lf->maxval);

    // Perform the histogram packing of the image
    for (v = 0; v < lf->v; v++) {
        for (u = 0; u < lf->u; u++) {
            for (t = 0; t < lf->t; t++) {
                for (s = 0; s < lf->s; s++) {
                    u_lf->val[v][u][t][s] = (img_t) backward_table[lf->val[v][u][t][s]];
                }
            }
        }
    }

    return u_lf;
}

// Write classes to file
void debug_predictors(DECODER *dec, int frame) {
    if (dec->debug_path != NULL) {
        char *predictor_file = (char *) alloc_mem(
                (strlen(dec->debug_path) + strlen("/predictors.txt\0") + 1) * sizeof(char));

        sprintf(predictor_file, "%s/predictors.txt", dec->debug_path);

        FILE *debug_fp = fopen(predictor_file, "a");

        fprintf(debug_fp, "SAI: %d\n", frame);
        fprintf(debug_fp, "#CLASS\t\tCOEFFICIENTS");

        for (int cl = 0; cl < dec->num_class; cl++) {
            fprintf(debug_fp, "\n%2d\t\t\t", cl);
            for (int p = 0; p < dec->full_prd_order; p++) {
                fprintf(debug_fp, " %3d", dec->predictor[cl][p]);
            }
        }

        fprintf(debug_fp, "\n");

        fclose(debug_fp);
        free(predictor_file);
    }
}

// Write partition to file
void debug_partition(DECODER *dec, int frame) {
    int d, k, i, j, z;
    unsigned short byte;

    if (dec->debug_path != NULL) {
        char *partition_img = (char *) alloc_mem(
                (strlen(dec->debug_path) + strlen("/partition_0000x0000_10bpp_LE_YUV444p.yuv") + 1) * sizeof(char));

        img_t ***qt_img = (img_t ***) alloc_3d_array(3, dec->ts[HEIGHT], dec->ts[WIDTH], sizeof(img_t));

        int elements = dec->ts[HEIGHT] * dec->ts[WIDTH] * 3 * (dec->depth == 8 ? 1 : 2);
        unsigned char *ptr, *qt_stream = (unsigned char *) alloc_mem(elements * sizeof(unsigned char));

        uint scale = (uint) floor(pow(2, dec->depth) / dec->num_class);

        for (i = 0; i < dec->ts[HEIGHT]; i++) {
            for (j = 0; j < dec->ts[WIDTH]; j++) {
                qt_img[0][i][j] = dec->org[0][i][j];

                qt_img[1][i][j] = dec->class[i][j] * scale;
                qt_img[2][i][j] = dec->class[i][j] * scale;
            }
        }

        if (dec->quadtree_depth > 0) {
            int blksize = MAX_BSIZE;

            // TODO: Rever e fazer bonito
            for (d = dec->quadtree_depth - 1; d >= 0; d--) {
                for (i = 0; i < dec->ts[HEIGHT]; i += blksize) {
                    for (j = 0; j < dec->ts[WIDTH]; j += blksize) {
                        for (k = 0; k < 2; k++) {
                            if ((dec->qtmap[d][i / blksize][j / blksize] == 0 && d == dec->quadtree_depth - 1) ||
                                (dec->qtmap[d][i / blksize][j / blksize] == 0 && d < dec->quadtree_depth - 1 &&
                                 dec->qtmap[d + 1][i / (blksize * 2)][j / (blksize * 2)] == 1)) {
                                if (i + blksize - 1 < dec->ts[HEIGHT]) {
                                    for (z = j; z < (j + blksize < dec->ts[WIDTH] ? j + blksize : dec->ts[WIDTH]); z++) {
                                        qt_img[k][i + blksize - 1][z] = (img_t) dec->maxval;
                                    }
                                }
                                if (j + blksize - 1 < dec->ts[WIDTH]) {
                                    for (z = i; z < (i + blksize < dec->ts[HEIGHT] ? i + blksize : dec->ts[HEIGHT]); z++) {
                                        qt_img[k][z][j + blksize - 1] = (img_t) dec->maxval;
                                    }
                                }
                            }
                            if (d == 0 && dec->qtmap[d][i / blksize][j / blksize] == 1) {
                                if (i + blksize / 2 - 1 < dec->ts[HEIGHT]) {
                                    for (z = j; z < (j + blksize < dec->ts[WIDTH] ? j + blksize : dec->ts[WIDTH]); z++) {
                                        qt_img[k][i + blksize / 2 - 1][z] = (img_t) dec->maxval;
                                    }
                                }
                                if (i + blksize - 1 < dec->ts[HEIGHT]) {
                                    for (z = j; z < (j + blksize < dec->ts[WIDTH] ? j + blksize : dec->ts[WIDTH]); z++) {
                                        qt_img[k][i + blksize - 1][z] = (img_t) dec->maxval;
                                    }
                                }
                                if (j + blksize / 2 - 1 < dec->ts[WIDTH]) {
                                    for (z = i; z < (i + blksize < dec->ts[HEIGHT] ? i + blksize : dec->ts[HEIGHT]); z++) {
                                        qt_img[k][z][j + blksize / 2 - 1] = (img_t) dec->maxval;
                                    }
                                }
                                if (j + blksize - 1 < dec->ts[WIDTH]) {
                                    for (z = i; z < (i + blksize < dec->ts[HEIGHT] ? i + blksize : dec->ts[HEIGHT]); z++) {
                                        qt_img[k][z][j + blksize - 1] = (img_t) dec->maxval;
                                    }
                                }
                            }
                        }
                    }
                }

                blksize >>= 1u;
            }
        }

        ptr = qt_stream;
        for (k = 0; k < 3; k++) {
            for (i = 0; i < dec->ts[HEIGHT]; i++) {
                for (j = 0; j < dec->ts[WIDTH]; j++) {
                    if (dec->depth > 8) {
                        byte = reverse_endianness(qt_img[k][i][j], LITTLE_ENDIANNESS);

                        *ptr++ = (byte >> 8u) & 0x00FFu;
                        *ptr++ = byte & 0x00FFu;
                    }
                    else {
                        *ptr++ = qt_img[k][i][j];
                    }
                }
            }
        }

        sprintf(partition_img, "%s/partition_%dx%d_%dbpp_LE_YUV444p.yuv", dec->debug_path, dec->ts[WIDTH], dec->ts[HEIGHT],
                dec->depth);

        FILE *debug_fp = fopen(partition_img, "ab");

        fwrite(qt_stream, sizeof(unsigned char), elements, debug_fp);

        fclose(debug_fp);

        free(qt_stream);
        free(partition_img);
        free(qt_img);

        // Quadtree map
        if (dec->quadtree_depth > 0) {
            char qtmap[1000];
            sprintf(qtmap, "%s/qtmap.txt", dec->debug_path);

            debug_fp = fopen(qtmap, "a");

            uint x, y, xx, yy;

            yy = (dec->ts[HEIGHT] + MAX_BSIZE - 1) / MAX_BSIZE;
            xx = (dec->ts[WIDTH] + MAX_BSIZE - 1) / MAX_BSIZE;

            fprintf(debug_fp, "SAI: %d\n", frame);
            fprintf(debug_fp, "#LEVEL\t\tPARTITION\n");

            for (i = dec->quadtree_depth - 1; i >= 0; i--) {
                fprintf(debug_fp, "%d\t\t\t", i);
                for (y = 0; y < yy; y++) {
                    for (x = 0; x < xx; x++) {
                        fprintf(debug_fp, " %d", dec->qtmap[i][y][x]);
                    }
                }

                fprintf(debug_fp, "\n");

                yy <<= 1u;
                xx <<= 1u;
            }

            fclose(debug_fp);
        }
    }
}

int main(int argc, char **argv) {
    int i, r;
    int vu[2], ts[2], num_class, maxval, depth, num_comp, num_group, endianness = LITTLE_ENDIANNESS;
    int num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta, hist_bytes;
    int *backward_table = NULL;
    int prd_order = 0;
    int sai_prd_order = 0;
    int sai_references, sai_distance_threshold;
    int format = SAI;
    int no_hilevels = 0, no_refsai_candidates = 0, no_refsai = 0, curr_frame = 0, curr_hilevel = 0, next_hilevel = 0;
    int **hilevels = NULL, *reference_list = NULL;
    int no_ra_regions = 0;
    bool use_current = false, use_disp = false;
    int disp_blk_size = 0, max_disparity = 0;
    SAIDISTANCE **reference_candidates = NULL;
    int sai_coord[2] = {0, 0};
    char *format_name = NULL;
    int debug = 0;
    char *debug_path = NULL;

    LF4D *lf = NULL, *save_err = NULL;
    DECODER *dec = NULL;
    char *infile, *outfile;
    FILE *fp;

    cpu_time();
    setbuf(stdout, 0);
    infile = outfile = NULL;

    for (i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch (argv[i][1]) {
                case 'E':
                    endianness = (int) strtol(argv[++i], NULL, 10);

                    if (endianness != LITTLE_ENDIANNESS && endianness != BIG_ENDIANNESS) {
                        endianness = LITTLE_ENDIANNESS;
                    }

                    break;

                case 'd':
                    debug = 1;

                    break;

                case 'r':
                    format_name = argv[++i];

                    if (strcmp(format_name, "MIA") == 0) {
                        format = MIA;
                        printf("\tMIA; --> Not yet implemented;\n");
                        exit(-2);
                    }
                    else if (strcmp(format_name, "PVS") == 0) {
                        format = PVS;
                    }
                    else if (strcmp(format_name, "SAI") == 0) {
                        format = SAI;
                    }
                    else {
                        printf("Light field format not recognised: %s.\n", format_name);
                        printf("Supported formats:\n");
                        printf("\tMIA; --> Not yet implemented;\n"
                               "\tPVS;\n"
                               "\tSAI.\n");
                        exit(-2);
                    }

                    break;

                default:
                    fprintf(stderr, "Unknown option: %s!\n", argv[i]);
                    exit(1);
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
        printf(BANNER"\n", MRP_VERSION, MRP_VERSION_DATE);
        printf("usage: decmrp [options] infile outfile\n");
        printf("-E num      Endianness: little-endian = 0, big-endian = 1. Default: %s\n", "little-endian");
        printf("-d          Create extra debug output (coefficients, partitions, etc.)\n");
        printf("-r str      Light field file format [%s]. Supported formats:\n", "SAI");
        printf("                MIA; --> Not yet implemented\n");
        printf("                PVS;\n");
        printf("                SAI.\n");
        printf("infile:     Input file\n");
        printf("outfile:    Output file\n");
        exit(0);
    }

    if (access(outfile, F_OK) != -1) {
        if (remove(outfile) != 0) {
            printf("Error deleting file %s.\n", outfile);
            return -1;
        }
    }

    fp = fileopen(infile, "rb");

    // Read header
    read_header(fp, vu, ts, &maxval, &depth, &num_comp, &num_group, &prd_order, &sai_prd_order, &num_pmodel,
                &coef_precision, &pm_accuracy, &f_huffman, &quadtree_depth, &delta, &hist_bytes, &no_hilevels,
                &no_ra_regions, &use_current, &use_disp, &disp_blk_size, &max_disparity, &sai_references, &sai_distance_threshold);

    // Alloc memory for the hierarchical information
    hilevels = (int **) alloc_2d_array(no_hilevels, vu[HEIGHT] * vu[WIDTH] + 1, sizeof(int));
    for (i = 0; i < no_hilevels; i++) hilevels[i][0] = 0;
    reference_candidates = (SAIDISTANCE **) alloc_mem((vu[HEIGHT] * vu[WIDTH]) * sizeof(SAIDISTANCE *));
    for (r = 0; r < vu[HEIGHT] * vu[WIDTH]; r++) reference_candidates[r] = (SAIDISTANCE *) alloc_mem(sizeof(SAIDISTANCE));
    reference_list = (int *) alloc_mem(sai_references * sizeof(int));

    // Create directory to store debug information
    if (debug == 1) {
        // Get debug path name
        int len = strlen(infile);
        debug_path = (char *) alloc_mem((len + 1) * sizeof(char));
        strcpy(debug_path, infile);
        strcpy(&debug_path[len - 4], "_dec");

        // Check if the directory exists and in that case deletes it
        struct stat sb;
        if (stat(debug_path, &sb) == 0 && (uint) S_IFDIR & sb.st_mode) {
            char remove_command[1000];
            sprintf(remove_command, "rm -r %s", debug_path);

            if (system(remove_command) == -1) {
                fprintf(stderr, "The program was unable to remove the debug path.\n");
                exit(EXIT_FAILURE);
            }
        }

        // Creates debug path
        if (mkdir(debug_path, 0775) == -1) {
            printf("Debug folder could not be created, the program will exit now.\n");
            exit(EXIT_FAILURE);
        }
    }

    if (hist_bytes != 0) {
        backward_table = decode_lookuptable(fp, hist_bytes, (int) pow(2, depth));
    }

    // Allocate memory for the decoded image
    lf = alloc_lf4d(vu[HEIGHT], vu[WIDTH], ts[HEIGHT], ts[WIDTH], (int) (pow(2, depth) - 1));
    // Allocate memory to save the error structure
    save_err = alloc_lf4d(vu[HEIGHT], vu[WIDTH], ts[HEIGHT], ts[WIDTH], (int) (pow(2, depth) - 1));

    printf("\nMRP-Video Decoder\n\n");
    // Print file characteristics to screen
    printf("%s (%d x %d x %d x %d x %d) -> %s\n", infile, vu[HEIGHT], vu[WIDTH], ts[HEIGHT], ts[WIDTH], depth, outfile);

    // Print coding parameters to screen
    printf("P = %d, V = %d, A = %d, D = %d\n\n", coef_precision, num_pmodel, pm_accuracy, delta);
    if (backward_table != NULL) {
        printf("Histogram packing on\n\n");
    }
    // Print prediction parameters to screen
    printf("LF Prediction order: ");
    printf("%d | ", prd_order);
    for (r = 0; r < sai_references; r++) printf("%d ", sai_prd_order);
    printf("(4DLF dimensions: %d x %d x %d x %d)\n\n", vu[HEIGHT], vu[WIDTH], ts[HEIGHT], ts[WIDTH]);

    // Loop all SAIs
    for (int s = 0; s < vu[HEIGHT] * vu[WIDTH]; s++) {
        num_class = read_class(fp);
        read_hilevel(&next_hilevel, &curr_frame, fp);

        printf("Decoding SAI %d, Hi Level %d", curr_frame, next_hilevel);

        if (curr_hilevel != next_hilevel) {
            no_refsai_candidates = 0;

            // Loop previous hierarchical levels to get the reference SAI candidates
            for (int hh = 0; hh < next_hilevel; hh++) {
                for (int f = 1; f <= hilevels[hh][0]; f++) {
                    // Store the reference SAI number on the reference_candidates structure
                    reference_candidates[no_refsai_candidates]->sai = hilevels[hh][f];

                    // Convert the SAI number to coordinates in the angular dimensions array
                    frame2coordinates(reference_candidates[no_refsai_candidates++]->coordinates, hilevels[hh][f], vu[WIDTH]);
                }
            }
        }

        curr_hilevel = next_hilevel;

        hilevels[curr_hilevel][0]++;
        hilevels[curr_hilevel][hilevels[curr_hilevel][0]] = curr_frame;

        no_refsai = 0;

        // Get the coordinates of the current SAI
        frame2coordinates(sai_coord, curr_frame, vu[WIDTH]);

        if (curr_hilevel > 0 && no_refsai_candidates > 0) {
            int curr_ra_region = 0, ref_ra_region = 0;

            // Calculate the distance of the reference candidates to the current SAI
            for (i = 0; i < no_refsai_candidates; i++) {
                reference_candidates[i]->distance = sqrt(pow(reference_candidates[i]->coordinates[0] - sai_coord[0], 2) + pow(reference_candidates[i]->coordinates[1] - sai_coord[1], 2));
            }

            // Sort the reference candidates by growing distance
            qsort(reference_candidates, no_refsai_candidates, sizeof(SAIDISTANCE *), compare_distance);

            // Select up to sai_references provided the distance is lower than sai_distance_threshold
            int reference_limit = no_refsai_candidates > sai_references ? sai_references : no_refsai_candidates;

            for (i = 0; i < reference_limit; i++) {
                if (no_ra_regions > 0) {
                    curr_ra_region = get_random_access_region(no_ra_regions, vu, sai_coord[HEIGHT], sai_coord[WIDTH],NULL);
                    ref_ra_region = get_random_access_region(no_ra_regions, vu, reference_candidates[i]->coordinates[HEIGHT],
                                                             reference_candidates[i]->coordinates[WIDTH], sai_coord);
                }

                if (reference_candidates[i]->distance < sai_distance_threshold &&
                    (no_ra_regions == 0 || (no_ra_regions > 0 && curr_ra_region == ref_ra_region))) {
                    reference_list[no_refsai++] = reference_candidates[i]->sai;
                }
            }

            printf(" (Refs:");
            for (r = 0; r < no_refsai; r++) printf(" %d", reference_list[r]);
            printf(")");
        }
        printf(":\n");

        if (use_current) {
            // Store the reference SAI number on the reference_candidates structure
            reference_candidates[no_refsai_candidates]->sai = curr_frame;

            // Convert the SAI number to coordinates in the angular dimensions array
            frame2coordinates(reference_candidates[no_refsai_candidates++]->coordinates, curr_frame, vu[WIDTH]);
        }

        //// Start of the decoding process
        printf("    Process started");

        // Creates new DECODER structure
        dec = init_decoder(fp, ts, maxval, num_class, num_comp, num_group, prd_order, sai_prd_order, coef_precision,
                           f_huffman, quadtree_depth, num_pmodel, pm_accuracy, delta, depth, curr_frame, lf, save_err,
                           no_refsai, reference_list, debug_path, use_disp, disp_blk_size, max_disparity);

        decode_class(fp, dec);

        decode_predictor(fp, dec);

        if (debug == 1) {
            debug_predictors(dec, curr_frame);
        }

        decode_threshold(fp, dec);

        dec->pmodels = init_pmodels(dec->num_group, dec->num_pmodel, dec->pm_accuracy, dec->pm_idx, dec->sigma, dec->maxval + 1);

        decode_image(fp, dec);

        if (debug == 1) {
            debug_partition(dec, curr_frame);
        }

        printf(" --> Process completed\n");

        // Save the decoded image and the residuals into the 4D structure
        for (int t = 0; t < dec->ts[HEIGHT]; t++) {
            for (int ss = 0; ss < dec->ts[WIDTH]; ss++) {
                lf->val[sai_coord[0]][sai_coord[1]][t][ss] = dec->org[0][t][ss];
                save_err->val[sai_coord[0]][sai_coord[1]][t][ss] = dec->err[0][t][ss];
            }
        }

        free_decoder(dec);
    }

    if (backward_table != NULL) {
        LF4D *unpacked = histogram_unpacking(lf, backward_table);

        write_yuv(unpacked, outfile, depth, endianness, format);

        safefree_lf4d(&unpacked);
    }
    else {
        write_yuv(lf, outfile, depth, endianness, format);
    }

    if (debug == 1) {
        free(debug_path);
    }

    for (r = 0; r < vu[HEIGHT] * vu[WIDTH]; r++) {
        free(reference_candidates[r]);
    }
    free(reference_candidates);
    free(reference_list);
    free(hilevels);

    safefree_lf4d(&lf);
    safefree_lf4d(&save_err);

    safefree((void **) &backward_table);

    fclose(fp);

    printf("\ncpu time: %.2f sec.\n", cpu_time());

    return(0);
}
