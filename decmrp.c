/*
“Commons Clause” License Condition v1.0

The Software is provided to you by the Licensor under the License, as defined below, subject to the following condition.

Without limiting other conditions in the License, the grant of rights under the License will not include, and the License does not grant to you, the right to Sell the Software.

For purposes of the foregoing, “Sell” means practicing any or all of the rights granted to you under the License to provide to third parties, for a fee or other consideration (including without limitation fees for hosting or consulting/ support services related to the Software), a product or service whose value derives, entirely or substantially, from the functionality of the Software. Any license notice or attribution required by the License must also include this Commons Clause License Condition notice.

Software: Dual Tree 4D Minimum Rate Predictors

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
#include <sys/stat.h>
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

DECODER *init_decoder(FILE *fp, const int *vu, const int *ts, int maxval, int num_comp, int num_group, int prd_order,
                      const int *mi_prd_order, int coef_precision, int f_huffman, int quadtree_depth, int num_pmodel,
                      int pm_accuracy, int delta, int depth, char *debug_path) {
    DECODER *dec;
    int i;

    dec = (DECODER *) alloc_mem(sizeof(DECODER));

    dec->vu[HEIGHT] = vu[HEIGHT];
    dec->vu[WIDTH] = vu[WIDTH];
    dec->ts[HEIGHT] = ts[HEIGHT];
    dec->ts[WIDTH] = ts[WIDTH];

    dec->maxval = maxval;
    dec->num_comp = num_comp;
    dec->num_group = num_group;

    dec->prd_order = prd_order;
    dec->mi_prd_order[UP] = mi_prd_order[UP];
    dec->mi_prd_order[LEFT] = mi_prd_order[LEFT];
    dec->mi_prd_order[RDIAG] = mi_prd_order[RDIAG];
    dec->mi_prd_order[LDIAG] = mi_prd_order[LDIAG];

    dec->full_prd_order = dec->prd_order + dec->mi_prd_order[UP] + dec->mi_prd_order[LEFT] + dec->mi_prd_order[RDIAG] +
                          dec->mi_prd_order[LDIAG];

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

    dec->num_class = read_class(fp);

    dec->predictor = (int **) alloc_2d_array(dec->num_class, dec->full_prd_order, sizeof(int));

    dec->err = (int ****) alloc_4d_array(dec->vu[HEIGHT] + 1, dec->vu[WIDTH], dec->ts[HEIGHT], dec->ts[WIDTH], sizeof(int));
    dec->err[dec->vu[HEIGHT]][0][0][0] = (dec->maxval + 1) >> 2;

    //Initiation of the reference offset
    dec->intra_roff = init_ref_offset(dec->vu, dec->ts, INTRA_PRED, dec->prd_order);

    // Keeps the weights used for the residue encoding context
    dec->intra_ctx_weight = init_ctx_weight(INTRA_PRED, dec->prd_order, dec->delta);

    if (dec->mi_prd_order[UP] > 0) {
        dec->mi_up_roff = init_ref_offset(dec->vu, dec->ts, MI_UP_PRED, dec->mi_prd_order[UP]);
        dec->mi_up_ctx_weight = init_ctx_weight(MI_UP_PRED, dec->mi_prd_order[UP], dec->delta);

        dec->mi_borders_ctx_weight = init_ctx_weight(MI_BORDERS_PRED, dec->mi_prd_order[UP], dec->delta);
    }

    if (dec->mi_prd_order[LEFT] > 0) {
        dec->mi_left_roff = init_ref_offset(dec->vu, dec->ts, MI_LEFT_PRED, dec->mi_prd_order[LEFT]);
        dec->mi_left_ctx_weight = init_ctx_weight(MI_LEFT_PRED, dec->mi_prd_order[LEFT], dec->delta);

        if (dec->mi_prd_order[LEFT] > dec->mi_prd_order[UP]) {
            dec->mi_borders_ctx_weight = init_ctx_weight(MI_BORDERS_PRED, dec->mi_prd_order[LEFT], dec->delta);
        }
    }

    if (dec->mi_prd_order[LDIAG] > 0) {
        dec->mi_ldiag_roff = init_ref_offset(dec->vu, dec->ts, MI_LDIAG_PRED, dec->mi_prd_order[LDIAG]);
        dec->mi_ldiag_ctx_weight = init_ctx_weight(MI_LDIAG_PRED, dec->mi_prd_order[LDIAG], dec->delta);

        if (dec->mi_prd_order[LDIAG] > dec->mi_prd_order[UP] && dec->mi_prd_order[LDIAG] > dec->mi_prd_order[LEFT]) {
            dec->mi_borders_ctx_weight = init_ctx_weight(MI_BORDERS_PRED, dec->mi_prd_order[LDIAG], dec->delta);
        }
    }

    if (dec->mi_prd_order[RDIAG] > 0) {
        dec->mi_rdiag_roff = init_ref_offset(dec->vu, dec->ts, MI_RDIAG_PRED, dec->mi_prd_order[RDIAG]);
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

    // Class and uquant arrays
    dec->class = (char ****) alloc_4d_array(dec->vu[HEIGHT], dec->vu[WIDTH], dec->ts[HEIGHT], dec->ts[WIDTH], sizeof(char));

    dec->uquant = (char **) alloc_2d_array(dec->num_class, MAX_UPARA + 1, sizeof(char));

    if (dec->num_pmodel > 1) {
        dec->pm_idx = (int *) alloc_mem(dec->num_group * sizeof(int));
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

    dec->mtfbuf = (int *) alloc_mem(dec->num_class * sizeof(int));

    return (dec);
}

void free_decoder(DECODER *dec) {
    int i, j, gr, num_subpm;

    free(dec->predictor);
    free(dec->err);

    free(dec->intra_ctx_weight);
    free_ref_offset(dec->vu, dec->ts, INTRA_PRED, dec->prd_order, dec->intra_roff);

    if(dec->mi_prd_order[UP] > 0) {
        free(dec->mi_up_ctx_weight);
        free_ref_offset(dec->vu, dec->ts, MI_UP_PRED, dec->mi_prd_order[UP], dec->mi_up_roff);
    }

    if(dec->mi_prd_order[LEFT] > 0) {
        free(dec->mi_left_ctx_weight);
        free_ref_offset(dec->vu, dec->ts, MI_LEFT_PRED, dec->mi_prd_order[LEFT], dec->mi_left_roff);
    }

    if(dec->mi_prd_order[LDIAG] > 0) {
        free(dec->mi_ldiag_ctx_weight);
        free_ref_offset(dec->vu, dec->ts, MI_LDIAG_PRED, dec->mi_prd_order[LDIAG], dec->mi_ldiag_roff);
    }

    if(dec->mi_prd_order[RDIAG] > 0) {
        free(dec->mi_rdiag_ctx_weight);
        free_ref_offset(dec->vu, dec->ts, MI_RDIAG_PRED, dec->mi_prd_order[RDIAG], dec->mi_rdiag_roff);
    }

    if (dec->mi_prd_order[UP] > 0 || dec->mi_prd_order[LEFT] > 0 || dec->mi_prd_order[LDIAG] > 0 ||
        dec->mi_prd_order[RDIAG] > 0)
        free(dec->mi_borders_ctx_weight);

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

    if (dec->debug_path != NULL) {
        free(dec->debug_path);
    }

    free(dec);
}

void read_header(FILE *fp, int vu[2], int ts[2], int *maxval, int *depth, int *num_comp, int *num_group, int *prd_order,
                 int mi_prd_order[4], int *num_pmodel, int *coef_precision, int *pm_accuracy, int *f_huffman,
                 int *quadtree_depth, int *delta, int *hist_bytes) {
    int i = 0;

    if (getbits(fp, 16) != MAGIC_NUMBER) {
        fprintf(stderr, "Not a compressed file!\n");
        exit(1);
    }

    vu[HEIGHT] = getbits(fp, 16);
    vu[WIDTH] = getbits(fp, 16);
    ts[HEIGHT] = getbits(fp, 16);
    ts[WIDTH] = getbits(fp, 16);
    *maxval = getbits(fp, 16);
    *depth = getbits(fp, 6);
    *num_comp = getbits(fp, 4);
    *num_group = getbits(fp, 6);

    *prd_order = getbits(fp, 8);

    for (i = 0; i < 4; i++) {
        mi_prd_order[i] = getbits(fp, 8);
    }

    *delta = getbits(fp, 8);
    *num_pmodel = getbits(fp, 6) + 1;
    *coef_precision = getbits(fp, 4) + 1;
    *pm_accuracy = getbits(fp, 4) - 1;
    *f_huffman = getbits(fp, 1);
    *quadtree_depth = (getbits(fp, 1)) ? QUADTREE_DEPTH : -1;
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

void decode_qtindex(FILE *fp, DECODER *dec, VLC *vlc, PMODEL *cpm, int tlv, int tlu, int tlt, int tls, int vu_blksize[2],
                    int ts_blksize, int width, int ts_level) {
    int i, cl, v, u, t, s, brv, bru, brt, brs, vu_blksize_part[4][2];
    PMODEL *pm;

    brv = (tlv + vu_blksize[HEIGHT] < dec->vu[HEIGHT]) ? (tlv + vu_blksize[HEIGHT]) : dec->vu[HEIGHT];
    bru = (tlu + vu_blksize[WIDTH] < dec->vu[WIDTH]) ? (tlu + vu_blksize[WIDTH]) : dec->vu[WIDTH];
    brt = (tlt + ts_blksize < dec->ts[HEIGHT]) ? (tlt + ts_blksize) : dec->ts[HEIGHT];
    brs = (tls + ts_blksize < dec->ts[WIDTH]) ? (tls + ts_blksize) : dec->ts[WIDTH];

    if (tlv >= brv || tlu >= bru || tlt >= brt || tls >= brs) return;

    if (vu_blksize[HEIGHT] > MIN_BSIZE_VU || vu_blksize[WIDTH] > MIN_BSIZE_VU || ts_level > 0) {
        if (dec->f_huffman == 1) {
            i = getbits(fp, 1);
        }
        else {
            pm = &dec->spm;
            i = rc_decode(fp, dec->rc, pm, 0, pm->size);
        }

        switch (i) {
            case 1:
                ts_blksize >>= 1;

                // v and u
                decode_qtindex(fp, dec, vlc, cpm, tlv, tlu, tlt, tls, vu_blksize, ts_blksize, width, ts_level - 1);
                decode_qtindex(fp, dec, vlc, cpm, tlv, tlu, tlt, tls + ts_blksize, vu_blksize, ts_blksize, width, ts_level - 1);
                decode_qtindex(fp, dec, vlc, cpm, tlv, tlu, tlt + ts_blksize, tls, vu_blksize, ts_blksize, width, ts_level - 1);
                decode_qtindex(fp, dec, vlc, cpm, tlv, tlu, tlt + ts_blksize, tls + ts_blksize, vu_blksize, ts_blksize, brs, ts_level - 1);

                return;
                break;

            case 2:
                if (vu_blksize[HEIGHT] > MIN_BSIZE_VU && vu_blksize[WIDTH] > MIN_BSIZE_VU) {
                    // v and u partition
                    vu_blksize_part[0][HEIGHT] = vu_blksize[HEIGHT] >> 1;
                    vu_blksize_part[0][WIDTH] = vu_blksize[WIDTH] >> 1;

                    vu_blksize_part[3][HEIGHT] = vu_blksize[HEIGHT] - vu_blksize_part[0][HEIGHT];
                    vu_blksize_part[3][WIDTH] = vu_blksize[WIDTH] - vu_blksize_part[0][WIDTH];

                    vu_blksize_part[1][HEIGHT] = vu_blksize_part[0][HEIGHT];
                    vu_blksize_part[1][WIDTH] = vu_blksize_part[3][WIDTH];

                    vu_blksize_part[2][HEIGHT] = vu_blksize_part[3][HEIGHT];
                    vu_blksize_part[2][WIDTH] = vu_blksize_part[0][WIDTH];

                    decode_qtindex(fp, dec, vlc, cpm, tlv, tlu, tlt, tls, vu_blksize_part[0], ts_blksize, width, ts_level);
                    decode_qtindex(fp, dec, vlc, cpm, tlv, tlu + vu_blksize_part[0][WIDTH], tlt, tls, vu_blksize_part[1], ts_blksize, width, ts_level);
                    decode_qtindex(fp, dec, vlc, cpm, tlv + vu_blksize_part[0][HEIGHT], tlu, tlt, tls, vu_blksize_part[2], ts_blksize, width, ts_level);
                    decode_qtindex(fp, dec, vlc, cpm, tlv + vu_blksize_part[0][HEIGHT], tlu + vu_blksize_part[0][WIDTH], tlt, tls, vu_blksize_part[3], ts_blksize, brs, ts_level);
                }
                else if (vu_blksize[HEIGHT] > MIN_BSIZE_VU && vu_blksize[WIDTH] <= MIN_BSIZE_VU) {
                    // v and u partition
                    vu_blksize_part[0][HEIGHT] = vu_blksize[HEIGHT] >> 1;
                    vu_blksize_part[0][WIDTH] = vu_blksize[WIDTH];

                    vu_blksize_part[3][HEIGHT] = vu_blksize[HEIGHT] - vu_blksize_part[0][HEIGHT];
                    vu_blksize_part[3][WIDTH] = vu_blksize[WIDTH];

                    decode_qtindex(fp, dec, vlc, cpm, tlv, tlu, tlt, tls, vu_blksize_part[0], ts_blksize, width, ts_level);
                    decode_qtindex(fp, dec, vlc, cpm, tlv + vu_blksize_part[0][HEIGHT], tlu, tlt, tls, vu_blksize_part[3], ts_blksize, brs, ts_level);
                }
                else if (vu_blksize[HEIGHT] <= MIN_BSIZE_VU && vu_blksize[WIDTH] > MIN_BSIZE_VU) {
                    // v and u partition
                    vu_blksize_part[0][HEIGHT] = vu_blksize[HEIGHT];
                    vu_blksize_part[0][WIDTH] = vu_blksize[WIDTH] >> 1;

                    vu_blksize_part[3][HEIGHT] = vu_blksize[HEIGHT];
                    vu_blksize_part[3][WIDTH] = vu_blksize[WIDTH] - vu_blksize_part[0][WIDTH];

                    decode_qtindex(fp, dec, vlc, cpm, tlv, tlu, tlt, tls, vu_blksize_part[0], ts_blksize, width, ts_level);
                    decode_qtindex(fp, dec, vlc, cpm, tlv, tlu + vu_blksize_part[0][WIDTH], tlt, tls, vu_blksize_part[3], ts_blksize, brs, ts_level);
                }

                return;
                break;

            default:
                break;
        }
    }

    i = rc_decode(fp, dec->rc, cpm, 0, cpm->size);
    mtf_classlabel(dec->class, dec->mtfbuf, tlv, tlu, tlt, tls, ts_blksize, width, dec->num_class);

    for (cl = 0; cl < dec->num_class; cl++) {
        if (dec->mtfbuf[cl] == i) break;
    }

    for (v = tlv; v < brv; v++) {
        for (u = tlu; u < bru; u++) {
            for (t = tlt; t < brt; t++) {
                for (s = tls; s < brs; s++) {
                    dec->class[v][u][t][s] = (char) cl;
                }
            }
        }
    }
}

void decode_class(FILE *fp, DECODER *dec) {
    int i, j, v, u, t, s, vu_blksize[2], ts_blksize, ts_level;
    VLC *vlc;
    PMODEL *pm, cpm[1];

    if (dec->quadtree_depth >= 0) {
        ts_level = dec->quadtree_depth;

        vu_blksize[HEIGHT] = MAX_BSIZE > dec->vu[HEIGHT] ? dec->vu[HEIGHT] : MAX_BSIZE;
        vu_blksize[WIDTH] = MAX_BSIZE > dec->vu[WIDTH] ? dec->vu[WIDTH] : MAX_BSIZE;

        ts_blksize = MAX_BSIZE; // 32
    }
    else {
        ts_level = 0;

        vu_blksize[HEIGHT] = BASE_BSIZE > dec->vu[HEIGHT] ? dec->vu[HEIGHT] : BASE_BSIZE;
        vu_blksize[WIDTH] = BASE_BSIZE > dec->vu[WIDTH] ? dec->vu[WIDTH] : BASE_BSIZE;

        ts_blksize = BASE_BSIZE; // 8
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
        double p, pp;
        int ctx;
        int qtree_code[2], mtf_code[MAX_CLASS];

        pm = &dec->spm;
        if (dec->quadtree_depth > 0) {
            set_spmodel(pm, 7, -1);

            for (ctx = 0; ctx < 2; ctx++) {
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
        if (vu_blksize[HEIGHT] > MIN_BSIZE_VU || vu_blksize[WIDTH] > MIN_BSIZE_VU || ts_level > 0) {
            pm->size = 3;

            p = qtree_prob[qtree_code[0]];
            pm->freq[0] = (uint) (p * (1 << 10));
            pp = qtree_prob[qtree_code[1]] * (1 - p);
            pm->freq[1] = (uint) (pp * (1 << 10));
            pp = (1 - qtree_prob[qtree_code[1]]) * (1 - p);
            pm->freq[2] = (uint) (pp * (1 << 10));

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

    for (v = 0; v < dec->vu[HEIGHT]; v += vu_blksize[HEIGHT]) {
        for (u = 0; u < dec->vu[WIDTH]; u += vu_blksize[WIDTH]) {
            for (t = 0; t < dec->ts[HEIGHT]; t += ts_blksize) {
                for (s = 0; s < dec->ts[WIDTH]; s += ts_blksize) {
                    decode_qtindex(fp, dec, vlc, cpm, v, u, t, s, vu_blksize, ts_blksize, dec->ts[WIDTH], ts_level);
                }
            }
        }
    }

    if (dec->f_huffman == 1) {
        free_vlc(vlc);
    }
    else {
        free(cpm->freq);
    }
}

int calc_udec(DECODER *dec, int v, int u, int t, int s) {
    int uu, k, *err_p;
    int *intra_wt_p;
    int *mi_up_wt_p = NULL, *mi_left_wt_p = NULL, *mi_ldiag_wt_p = NULL, *mi_rdiag_wt_p = NULL;
    int *intra_roff_p = NULL;
    int *mi_up_roff_p = NULL, *mi_left_roff_p = NULL, *mi_ldiag_roff_p = NULL, *mi_rdiag_roff_p = NULL;

    err_p = &dec->err[v][u][t][s];

    intra_roff_p = dec->intra_roff[v][u][t][s];
    intra_wt_p = dec->intra_ctx_weight;

    if (dec->mi_prd_order[UP] > 0) {
        mi_up_roff_p = dec->mi_up_roff[v][u][t][s];

        mi_up_wt_p = (v < 1) ? dec->mi_borders_ctx_weight : dec->mi_up_ctx_weight;
    }

    if (dec->mi_prd_order[LEFT] > 0) {
        mi_left_roff_p = dec->mi_left_roff[v][u][t][s];

        mi_left_wt_p = (u < 1) ? dec->mi_borders_ctx_weight : dec->mi_left_ctx_weight;
    }

    if (dec->mi_prd_order[LDIAG] > 0) {
        mi_ldiag_roff_p = dec->mi_ldiag_roff[v][u][t][s];

        mi_ldiag_wt_p = (v < 1 || u < 1) ? dec->mi_borders_ctx_weight : dec->mi_ldiag_ctx_weight;
    }

    if (dec->mi_prd_order[RDIAG] > 0) {
        mi_rdiag_roff_p = dec->mi_rdiag_roff[v][u][t][s];
        mi_rdiag_wt_p = (v < 1 || u > dec->vu[WIDTH] - 1) ? dec->mi_borders_ctx_weight : dec->mi_rdiag_ctx_weight;
    }

    uu = 0;

    for (k = 0; k < dec->prd_order; k++) {
        uu += err_p[*intra_roff_p++] * (*intra_wt_p++);
    }

    for (k = 0; k < dec->mi_prd_order[UP]; k++) {
        uu += err_p[*mi_up_roff_p++] * (*mi_up_wt_p++);
    }

    for (k = 0; k < dec->mi_prd_order[LEFT]; k++) {
        uu += err_p[*mi_left_roff_p++] * (*mi_left_wt_p++);
    }

    for (k = 0; k < dec->mi_prd_order[LDIAG]; k++) {
        uu += err_p[*mi_ldiag_roff_p++] * (*mi_ldiag_wt_p++);
    }

    for (k = 0; k < dec->mi_prd_order[RDIAG]; k++) {
        uu += err_p[*mi_rdiag_roff_p++] * (*mi_rdiag_wt_p++);
    }

    uu >>= 6;

    if (uu > MAX_UPARA) uu = MAX_UPARA;

    return (uu);
}

int calc_prd(img_t ****org, DECODER *dec, int cl, int v, int u, int t, int s) {
    int k, prd, *coef_p;
    img_t *org_p;
    int *intra_roff_p = NULL;
    int *mi_up_roff_p = NULL, *mi_left_roff_p = NULL, *mi_ldiag_roff_p = NULL, *mi_rdiag_roff_p = NULL;

    intra_roff_p = dec->intra_roff[v][u][t][s];
    if (dec->mi_prd_order[UP] > 0) mi_up_roff_p = dec->mi_up_roff[v][u][t][s];
    if (dec->mi_prd_order[LEFT] > 0) mi_left_roff_p = dec->mi_left_roff[v][u][t][s];
    if (dec->mi_prd_order[LDIAG] > 0) mi_ldiag_roff_p = dec->mi_ldiag_roff[v][u][t][s];
    if (dec->mi_prd_order[RDIAG] > 0) mi_rdiag_roff_p = dec->mi_rdiag_roff[v][u][t][s];

    org_p = &org[v][u][t][s];
    coef_p = dec->predictor[cl];

    prd = 0;

    for (k = 0; k < dec->prd_order; k++) {
        prd += org_p[*intra_roff_p++] * (*coef_p++);
    }

    for (k = 0; k < dec->mi_prd_order[UP]; k++) {
        prd += org_p[*mi_up_roff_p++] * (*coef_p++);
    }

    for (k = 0; k < dec->mi_prd_order[LEFT]; k++) {
        prd += org_p[*mi_left_roff_p++] * (*coef_p++);
    }

    for (k = 0; k < dec->mi_prd_order[LDIAG]; k++) {
        prd += org_p[*mi_ldiag_roff_p++] * (*coef_p++);
    }

    for (k = 0; k < dec->mi_prd_order[RDIAG]; k++) {
        prd += org_p[*mi_rdiag_roff_p++] * (*coef_p++);
    }

    if (prd < 0) prd = 0;
    else if (prd > dec->maxprd) prd = dec->maxprd;

    return (prd);
}

LF4D *decode_image(FILE *fp, DECODER *dec) {
    int v, u, t, s, cl, gr, prd, uu, e, E, p;
    LF4D *lf = NULL;

    // TODO: consider add org to decoder
    img_t ****org = (img_t ****) alloc_4d_array(dec->vu[HEIGHT] + 1, dec->vu[WIDTH], dec->ts[HEIGHT], dec->ts[WIDTH],
                                                sizeof(img_t));
    org[dec->vu[HEIGHT]][0][0][0] = (img_t) (dec->maxval + 1) >> 1;

    if (dec->f_huffman == 1) {
        VLC *vlc;
        dec->vlcs = init_vlcs(dec->pmodels, dec->num_group, 1);

        for (v = 0; v < dec->vu[HEIGHT]; v++) {
            for (u = 0; u < dec->vu[WIDTH]; u++) {
                for (t = 0; t < dec->ts[HEIGHT]; t++) {
                    for (s = 0; s < dec->ts[WIDTH]; s++) {
                        cl = dec->class[v][u][t][s];
                        uu = calc_udec(dec, v, u, t, s);
                        gr = dec->uquant[cl][uu];
                        prd = calc_prd(org, dec, cl, v, u, t, s);
                        prd >>= (dec->coef_precision - 1);
                        p = (prd + 1) >> 1;
                        vlc = &dec->vlcs[gr][0];
                        dec->err[v][u][t][s] = E = decode_vlc(fp, vlc);
                        e = E2e(E, p, prd & 1, dec->maxval);
                        org[v][u][t][s] = (img_t) (p + e);
                    }
                }
            }
        }
    }
    else {
        PMODEL *pm;
        if (dec->pm_accuracy < 0) {
            for (v = 0; v < dec->vu[HEIGHT]; v++) {
                for (u = 0; u < dec->vu[WIDTH]; u++) {
                    for (t = 0; t < dec->ts[HEIGHT]; t++) {
                        for (s = 0; s < dec->ts[WIDTH]; s++) {
                            cl = dec->class[v][u][t][s];
                            uu = calc_udec(dec, v, u, t, s);
                            gr = dec->uquant[cl][uu];
                            prd = calc_prd(org, dec, cl, v, u, t, s);
                            prd >>= (dec->coef_precision - 1);
                            p = (prd + 1) >> 1;
                            pm = dec->pmodels[gr][0];
                            dec->err[v][u][t][s] = E = rc_decode(fp, dec->rc, pm, 0, pm->size);
                            e = E2e(E, p, prd & 1, dec->maxval);
                            org[v][u][t][s] = (img_t) (p + e);
                        }
                    }
                }
            }
        }
        else {
            int mask, shift, base;
            mask = (1 << dec->pm_accuracy) - 1;
            shift = dec->coef_precision - dec->pm_accuracy;

            for (v = 0; v < dec->vu[HEIGHT]; v++) {
                for (u = 0; u < dec->vu[WIDTH]; u++) {
                    for (t = 0; t < dec->ts[HEIGHT]; t++) {
                        for (s = 0; s < dec->ts[WIDTH]; s++) {
                            cl = dec->class[v][u][t][s];
                            uu = calc_udec(dec, v, u, t, s);
                            gr = dec->uquant[cl][uu];
                            prd = calc_prd(org, dec, cl, v, u, t, s);
                            base = (dec->maxprd - prd + (1 << shift) / 2) >> shift;
                            pm = dec->pmodels[gr][0] + (base & mask);
                            base >>= dec->pm_accuracy;
                            p = rc_decode(fp, dec->rc, pm, base, base + dec->maxval + 1) - base;
                            org[v][u][t][s] = (img_t) p;
                            prd >>= (dec->coef_precision - 1);
                            e = (p << 1) - prd;
                            dec->err[v][u][t][s] = (e > 0) ? (e - 1) : (-e);
                        }
                    }
                }
            }
        }
    }

    lf = alloc_lf4d(dec->vu[HEIGHT], dec->vu[WIDTH], dec->ts[HEIGHT], dec->ts[WIDTH], dec->maxval);

    for (v = 0; v < dec->vu[HEIGHT]; v++) {
        for (u = 0; u < dec->vu[WIDTH]; u++) {
            for (t = 0; t < dec->ts[HEIGHT]; t++) {
                for (s = 0; s < dec->ts[WIDTH]; s++) {
                    lf->val[v][u][t][s] = org[v][u][t][s];
                }
            }
        }
    }

    free(org);

    return (lf);
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
void debug_predictors(DECODER *dec) {
    if (dec->debug_path != NULL) {
        char *predictor_file = (char *) alloc_mem(
                (strlen(dec->debug_path) + strlen("/predictors.txt\0") + 1) * sizeof(char));

        sprintf(predictor_file, "%s/predictors.txt", dec->debug_path);

        FILE *debug_fp = fopen(predictor_file, "w");

        fprintf(debug_fp, "#CLASS\t\tCOEFFICIENTS");

        for (int cl = 0; cl < dec->num_class; cl++) {
            fprintf(debug_fp, "\n%2d\t\t\t", cl);
            for (int p = 0; p < dec->full_prd_order; p++) {
                fprintf(debug_fp, " %3d", dec->predictor[cl][p]);
            }
        }

        fclose(debug_fp);
        free(predictor_file);
    }
}
// Write classes to file
void debug_class(DECODER *dec) {
    if (dec->debug_path != NULL) {
        char *classes_file = (char *) alloc_mem(
                (strlen(dec->debug_path) + strlen("/classes.txt\0") + 1) * sizeof(char));

        sprintf(classes_file, "%s/classes.txt", dec->debug_path);

        FILE *debug_fp = fopen(classes_file, "w");

        fprintf(debug_fp, "#Pixel\t\t\tClass");

        for (int v = 0; v < dec->vu[HEIGHT]; v++) {
            for (int u = 0; u < dec->vu[WIDTH]; u++) {
                for (int t = 0; t < dec->ts[HEIGHT]; t++) {
                    for (int s = 0; s < dec->ts[WIDTH]; s++) {
                        fprintf(debug_fp, "%02d,%02d,%03d,%03d:\t\t%3d\n", v, u, t, s, dec->class[v][u][t][s]);
                    }
                }
            }
        }

        fclose(debug_fp);
        free(classes_file);
    }
}

// Write partition to file
void debug_partition(DECODER *dec, LF4D *lf, int endianness) {
//    int d, i, j, k, l, v, u, t, s;
//
//    img_t *qt_lf_ptr[3];
//
//    if (dec->debug_path != NULL) {
//        char *partition_img = (char *) alloc_mem(
//                (strlen(dec->debug_path) + strlen("/partition_0000x0000_10bpp_LE_YUV444p.yuv\0")) * sizeof(char));
//
//        LF4D *qt_lf[3];
//        for (k = 0; k < 3; k++) {
//            qt_lf[k] = alloc_lf4d(dec->vu[HEIGHT], dec->vu[WIDTH], dec->ts[HEIGHT], dec->ts[WIDTH], (int) pow(2, dec->depth) - 1);
//            qt_lf_ptr[k] = &qt_lf[k]->val[0][0][0][0];
//        }
//
//        uint scale = (uint) floor(pow(2, dec->depth) / dec->num_class);
//
//        for (v = 0; v < dec->vu[HEIGHT]; v++) {
//            for (u = 0; u < dec->vu[WIDTH]; u++) {
//                for (t = 0; t < dec->ts[HEIGHT]; t++) {
//                    for (s = 0; s < dec->ts[WIDTH]; s++) {
//                        *qt_lf_ptr[0]++ = lf->val[v][u][t][s];
//
//                        *qt_lf_ptr[1]++ = dec->class[v][u][t][s] * scale;
//                        *qt_lf_ptr[2]++ = dec->class[v][u][t][s] * scale;
//                    }
//                }
//            }
//        }
//
//        if (dec->quadtree_depth > 0) {
//            uint blksize = MAX_BSIZE;
//
//            // TODO: Rever e tentar fazer bonito
//            for (d = dec->quadtree_depth - 1; d >= 0; d--) {
//                for (v = 0; v < dec->vu[HEIGHT]; v += blksize) {
//                    for (u = 0; u < dec->vu[WIDTH]; u += blksize) {
//                        for (t = 0; t < dec->ts[HEIGHT]; t += blksize) {
//                            for (s = 0; s < dec->ts[WIDTH]; s += blksize) {
//                                for (k = v; k < (v + blksize < dec->vu[HEIGHT] ? v + blksize : dec->vu[HEIGHT]); k++) {
//                                    for (l = u;
//                                         l < (u + blksize < dec->vu[WIDTH] ? u + blksize : dec->vu[WIDTH]); l++) {
//                                        if ((dec->qtmap[d][v / blksize][u / blksize][t / blksize][s / blksize] == 0 &&
//                                             d == dec->quadtree_depth - 1) ||
//                                            (dec->qtmap[d][v / blksize][u / blksize][t / blksize][s / blksize] == 0 &&
//                                             d < dec->quadtree_depth - 1 &&
//                                             dec->qtmap[d + 1][v / (blksize * 2)][u / (blksize * 2)][t / (blksize * 2)][
//                                                     s / (blksize * 2)] == 1)) {
//                                            if (t + blksize - 1 < dec->ts[HEIGHT]) {
//                                                for (j = s; j < (s + blksize < dec->ts[WIDTH] ? s + blksize
//                                                                                              : dec->ts[WIDTH]); j++) {
//                                                    qt_lf[0]->val[k][l][t + blksize - 1][j] = (img_t) qt_lf[0]->maxval;
//                                                    qt_lf[1]->val[k][l][t + blksize - 1][j] = (img_t) qt_lf[0]->maxval;
//                                                    qt_lf[2]->val[k][l][t + blksize - 1][j] = (img_t) qt_lf[0]->maxval;
//                                                }
//                                            }
//                                            if (s + blksize - 1 < dec->ts[WIDTH]) {
//                                                for (i = t; i < (t + blksize < dec->ts[HEIGHT] ? t + blksize
//                                                                                               : dec->ts[HEIGHT]); i++) {
//                                                    qt_lf[0]->val[k][l][i][s + blksize - 1] = (img_t) qt_lf[0]->maxval;
//                                                    qt_lf[1]->val[k][l][i][s + blksize - 1] = (img_t) qt_lf[0]->maxval;
//                                                    qt_lf[2]->val[k][l][i][s + blksize - 1] = (img_t) qt_lf[0]->maxval;
//                                                }
//                                            }
//                                        }
//                                        if (d == 0 &&
//                                            dec->qtmap[d][v / blksize][u / blksize][t / blksize][s / blksize] == 1) {
//                                            if (t + blksize / 2 - 1 < dec->ts[HEIGHT]) {
//                                                for (j = s; j < (s + blksize < dec->ts[WIDTH] ? s + blksize
//                                                                                              : dec->ts[WIDTH]); j++) {
//                                                    qt_lf[0]->val[k][l][t + blksize / 2 - 1][j] = (img_t) qt_lf[0]->maxval;
//                                                    qt_lf[1]->val[k][l][t + blksize / 2 - 1][j] = (img_t) qt_lf[0]->maxval;
//                                                    qt_lf[2]->val[k][l][t + blksize / 2 - 1][j] = (img_t) qt_lf[0]->maxval;
//                                                }
//                                            }
//                                            if (t + blksize - 1 < dec->ts[HEIGHT]) {
//                                                for (j = s; j < (s + blksize < dec->ts[WIDTH] ? s + blksize
//                                                                                              : dec->ts[WIDTH]); j++) {
//                                                    qt_lf[0]->val[k][l][t + blksize - 1][j] = (img_t) qt_lf[0]->maxval;
//                                                    qt_lf[1]->val[k][l][t + blksize - 1][j] = (img_t) qt_lf[0]->maxval;
//                                                    qt_lf[2]->val[k][l][t + blksize - 1][j] = (img_t) qt_lf[0]->maxval;
//                                                }
//                                            }
//                                            if (s + blksize / 2 - 1 < dec->ts[WIDTH]) {
//                                                for (i = t; i < (t + blksize < dec->ts[HEIGHT] ? t + blksize
//                                                                                               : dec->ts[HEIGHT]); i++) {
//                                                    qt_lf[0]->val[k][l][i][s + blksize / 2 - 1] = (img_t) qt_lf[0]->maxval;
//                                                    qt_lf[1]->val[k][l][i][s + blksize / 2 - 1] = (img_t) qt_lf[0]->maxval;
//                                                    qt_lf[2]->val[k][l][i][s + blksize / 2 - 1] = (img_t) qt_lf[0]->maxval;
//                                                }
//                                            }
//                                            if (s + blksize - 1 < dec->ts[WIDTH]) {
//                                                for (i = t; i < (t + blksize < dec->ts[HEIGHT] ? t + blksize
//                                                                                               : dec->ts[HEIGHT]); i++) {
//                                                    qt_lf[0]->val[k][l][i][s + blksize - 1] = (img_t) qt_lf[0]->maxval;
//                                                    qt_lf[1]->val[k][l][i][s + blksize - 1] = (img_t) qt_lf[0]->maxval;
//                                                    qt_lf[2]->val[k][l][i][s + blksize - 1] = (img_t) qt_lf[0]->maxval;
//                                                }
//                                            }
//                                        }
//                                    }
//                                }
//                            }
//                        }
//                    }
//                }
//
//                blksize >>= 1u;
//            }
//
//            for (v = 0; v < dec->vu[HEIGHT]; v++) {
//                for (u = 0; u < dec->vu[WIDTH]; u++) {
//                    for (t = 0; t < dec->ts[HEIGHT]; t++) {
//                        qt_lf[0]->val[v][u][t][dec->ts[WIDTH] - 1] = (img_t) qt_lf[0]->maxval;
//                        qt_lf[1]->val[v][u][t][dec->ts[WIDTH] - 1] = (img_t) qt_lf[0]->maxval;
//                        qt_lf[2]->val[v][u][t][dec->ts[WIDTH] - 1] = (img_t) qt_lf[0]->maxval;
//                    }
//                    for (s = 0; s < dec->ts[WIDTH]; s++) {
//                        qt_lf[0]->val[v][u][dec->ts[HEIGHT] - 1][s] = (img_t) qt_lf[0]->maxval;
//                        qt_lf[1]->val[v][u][dec->ts[HEIGHT] - 1][s] = (img_t) qt_lf[0]->maxval;
//                        qt_lf[2]->val[v][u][dec->ts[HEIGHT] - 1][s] = (img_t) qt_lf[0]->maxval;
//                    }
//                }
//            }
//        }
//
//        // Quadtree partition image name
//        sprintf(partition_img, "%s/partition_%dx%d_%dbpp_LE_YUV444p.yuv", dec->debug_path,
//                dec->vu[WIDTH] * dec->ts[WIDTH], dec->vu[HEIGHT] * dec->ts[HEIGHT], dec->depth);
//
//        // Write image to file
//        for (k = 0; k < 3; k++) {
//            write_yuv(qt_lf[k], partition_img, dec->depth, endianness, SAI);
//        }
//
//        free(partition_img);
//
//        for (k = 0; k < 3; k++) {
//            safefree_lf4d(&qt_lf[k]);
//        }
//    }
}

// Write residuals to file
void debug_err(DECODER *dec, LF4D *lf, int endianness) {
    if (dec->debug_path != NULL) {
        char *residuals_img = (char *) alloc_mem(
                (strlen(dec->debug_path) + strlen("/residuals_0000x0000_10bpp_LE_YUV444p.yuv\0")) * sizeof(char));

        // Quadtree partition image name
        sprintf(residuals_img, "%s/residuals_%dx%d_%dbpp_LE_GRAY.yuv", dec->debug_path,
                dec->vu[WIDTH] * dec->ts[WIDTH], dec->vu[HEIGHT] * dec->ts[HEIGHT], dec->depth + 1);

        LF4D *residuals;
        residuals = alloc_lf4d(dec->vu[HEIGHT], dec->vu[WIDTH], dec->ts[HEIGHT], dec->ts[WIDTH], (int) pow(2, dec->depth + 1) - 1);

        int v, u, t, s, cl, prd, e, p;
        int mask, shift, base;

        for (v = 0; v < dec->vu[HEIGHT]; v++) {
            for (u = 0; u < dec->vu[WIDTH]; u++) {
                for (t = 0; t < dec->ts[HEIGHT]; t++) {
                    for (s = 0; s < dec->ts[WIDTH]; s++) {
                        cl = dec->class[v][u][t][s];
                        prd = calc_prd(lf->val, dec, cl, v, u, t, s);
                        prd >>= (dec->coef_precision - 1);

                        residuals->val[v][u][t][s] = (((int) lf->val[v][u][t][s] << 1) - prd >> 1) + (int) pow(2, dec->depth) - 1;
                    }
                }
            }
        }

        // Write image to file
        write_yuv(residuals, residuals_img, dec->depth + 1, endianness, SAI);

        free(residuals_img);
        safefree_lf4d(&residuals);
    }
}

int main(int argc, char **argv) {
    int i;
    int vu[2], ts[2], maxval, depth, num_comp, num_group, endianness = LITTLE_ENDIANNESS;
    int num_pmodel, coef_precision, pm_accuracy, f_huffman, quadtree_depth, delta, hist_bytes;
    int *backward_table = NULL;
    int prd_order = 0;
    int mi_prd_order[4] = {0, 0, 0, 0};
    int format = SAI;
    char *format_name = NULL;
    int debug = 0;
    char *debug_path = NULL;

    LF4D *lf = NULL;
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
        printf("                MIA; --> Not yet implemented\n"
               "                PVS;\n"
               "                SAI.\n");
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
    read_header(fp, vu, ts, &maxval, &depth, &num_comp, &num_group, &prd_order, mi_prd_order, &num_pmodel,
                &coef_precision, &pm_accuracy, &f_huffman, &quadtree_depth, &delta, &hist_bytes);

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
            char aux[1000];
            sprintf(aux, "rm -r %s", debug_path);
            (void) system(aux);
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

    printf("\nMRP-Video Decoder\n\n");
    // Print file characteristics to screen
    printf("%s (%d x %d x %d x %d x %d) -> %s\n", infile, vu[HEIGHT], vu[WIDTH], ts[HEIGHT], ts[WIDTH], depth, outfile);

    // Print coding parameters to screen
    printf("P = %d, V = %d, A = %d, D = %d\n\n", coef_precision, num_pmodel, pm_accuracy, delta);
    if (backward_table != NULL) {
        printf("Histogram packing on\n\n");
    }
    // Print prediction parameters to screen
    printf("Prediction order:\n\tFrame I: %d\n\n", prd_order);

    printf("LF Prediction order: %d %d %d %d ", mi_prd_order[UP], mi_prd_order[LEFT], mi_prd_order[LDIAG],
           mi_prd_order[RDIAG]);
    printf("(4DLF dimensions: %d x %d x %d x %d)\n\n", vu[HEIGHT], vu[WIDTH], ts[HEIGHT], ts[WIDTH]);

    //// Start of the decoding process
    printf("Process started");

    // Creates new DECODER structure
    dec = init_decoder(fp, vu, ts, maxval, num_comp, num_group, prd_order, mi_prd_order, coef_precision, f_huffman,
                       quadtree_depth, num_pmodel, pm_accuracy, delta, depth, debug_path);

    decode_class(fp, dec);

    decode_predictor(fp, dec);

    if (debug == 1) {
        debug_predictors(dec);
        debug_class(dec);
    }

    decode_threshold(fp, dec);

    dec->pmodels = init_pmodels(dec->num_group, dec->num_pmodel, dec->pm_accuracy, dec->pm_idx, dec->sigma, dec->maxval + 1);

    lf = decode_image(fp, dec);

    if (debug == 1) {
        debug_partition(dec, lf, endianness);
        debug_err(dec, lf, endianness);

        free(debug_path);
    }

    printf(" --> Process completed\n");

    if (backward_table != NULL) {
        LF4D *unpacked = histogram_unpacking(lf, backward_table);

        write_yuv(unpacked, outfile, depth, endianness, format);

        safefree_lf4d(&unpacked);
    }
    else {
        write_yuv(lf, outfile, depth, endianness, format);
    }

    free_decoder(dec);
    safefree_lf4d(&lf);

    safefree((void **) &backward_table);

    fclose(fp);

    printf("\ncpu time: %.2f sec.\n", cpu_time());

    return(0);
}
