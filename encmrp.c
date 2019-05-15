#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <libgen.h>
#include <math.h>
#include <limits.h>
#include <sys/stat.h>
#include "mrp.h"
#include "mrp_config.h"

extern POINT dyx[];
extern POINT idyx[];
extern POINT tridyx[];
extern double sigma_h[], sigma_a[];
extern double qtree_prob[];

// TODO: fix functions headers
/*---------------------------- init_encoder ------------------------*
 |  Function init_encoder
 |
 |  Purpose:  Initializes the encoder structure
 |
 |  Parameters:
 |      img		 		--> Image structure with the frame to encode (IN)
 |      back_ref_img	--> Image structure with the backward reference frame (when available) (IN)
 |      for_ref_img	    --> Image structure with the forward reference frame (when available) (IN)
 |      back_ref_error	--> Array with the error of the backward reference frame (when available) (IN)
 |      for_ref_error   --> Array with the error of the forward reference frame (when available) (IN)
 |		num_class		--> Number of classes to use (number of different predictors) (IN)
 |		num_group		--> Number of pixels that are taken into account to form the context for entropy coding of the residue (IN)
 |		prd_order		--> Order of the predictors (number of pixels to use) (IN)
 |		back_prd_order	--> Order of the predictors in the backward reference frame (number of pixels to use, when available) (IN)
 |		for_prd_order	--> Order of the predictors in the forward reference frame (number of pixels to use, when available) (IN)
 |		mi_prd_order	--> Order of the predictors in neighbour micro images (number of pixels to use, when available) (IN)
 |		coef_precision	--> Precision of the coefficients (IN)
 |		f_huffman		--> Huffman coding flag (IN)
 |		quadtree_depth	--> Quadtree flag (IN)
 |		num_pmodel		--> Number of probability models (IN)
 |		pm_accuracy		--> Probability model accuracy (IN)
 |		delta			--> Parameter representing the distance between frames (IN)
 |		depth			--> Input sequence bit depth (IN)
 |
 |  Returns:  ENCODER* --> returns a encoder type structure
 *-------------------------------------------------------------------*/
ENCODER *
init_encoder(LF4D *lf, int num_class, int num_group, int prd_order, const int mi_prd_order[4], int coef_precision,
             int f_huffman, int quadtree_depth, int num_pmodel, int pm_accuracy, int delta, int depth,
             char *debug_path) {
    //Declare new encoder struct
    ENCODER *enc;
    int i, j, v, u, t, s;
    double c;

    //Allocation of the memory for the encoder struct
    enc = (ENCODER *) alloc_mem(sizeof(ENCODER));

    //Copy of the video/image properties to encoder
    enc->vu[HEIGHT] = lf->v;
    enc->vu[WIDTH] = lf->u;
    enc->ts[HEIGHT] = lf->t;
    enc->ts[WIDTH] = lf->s;
    enc->maxval = lf->maxval;

    // Copy of the encoding parameters
    enc->num_class = num_class; // M
    enc->num_group = num_group; // 16

    enc->prd_order = prd_order; // K
    enc->mi_prd_order[UP] = mi_prd_order[UP];
    enc->mi_prd_order[LEFT] = mi_prd_order[LEFT];
    enc->mi_prd_order[LDIAG] = mi_prd_order[LDIAG];
    enc->mi_prd_order[RDIAG] = mi_prd_order[RDIAG];

    enc->full_prd_order = enc->prd_order + enc->mi_prd_order[UP] + enc->mi_prd_order[LEFT] + enc->mi_prd_order[LDIAG] +
                          enc->mi_prd_order[RDIAG];

    enc->coef_precision = coef_precision; // P
    enc->f_huffman = f_huffman; // h
    enc->quadtree_depth = quadtree_depth; // f
    enc->num_pmodel = num_pmodel; // V
    enc->pm_accuracy = pm_accuracy; // A
    enc->maxprd = enc->maxval << enc->coef_precision; // Maximum prediction value allowed
    enc->delta = delta; // Parameter representing the distance between frames
    enc->depth = depth;
    enc->etype = 0;

    if (debug_path != NULL) {
        enc->debug_path = (char *) alloc_mem((strlen(debug_path) + 1) * sizeof(char));
        strcpy(enc->debug_path, debug_path);
    }
    else {
        enc->debug_path = NULL;
    }

    // Alloc memory to predictors array
    enc->predictor = (int **) alloc_2d_array(enc->num_class, enc->full_prd_order, sizeof(int));

    // Alloc memory to array of the threshold values used to quantize the weighted sum of neighboring residues
    enc->th = (int **) alloc_2d_array(enc->num_class, enc->num_group, sizeof(int));

    // enc->th initializing cycle
    for (i = 0; i < enc->num_class; i++) {
        for (j = 0; j < enc->num_group - 1; j++) {
            enc->th[i][j] = 0;
        }
        enc->th[i][enc->num_group - 1] = MAX_UPARA + 1;
    }

    // More memory allocation
    enc->upara = (int ****) alloc_4d_array(enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH], sizeof(int));
    enc->prd = (int ****) alloc_4d_array(enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH], sizeof(int));
    // Memory allocation for the original image
    enc->org = (int ****) alloc_4d_array(enc->vu[HEIGHT] + 1, enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH], sizeof(int));
    // Memory allocation for the error image
    enc->err = (int ****) alloc_4d_array(enc->vu[HEIGHT] + 1, enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH], sizeof(int));

    // Quadtree map
    if (enc->quadtree_depth > 0) {
        v = (enc->vu[HEIGHT] + MAX_BSIZE - 1) / MAX_BSIZE;
        u = (enc->vu[WIDTH] + MAX_BSIZE - 1) / MAX_BSIZE;
        t = (enc->ts[HEIGHT] + MAX_BSIZE - 1) / MAX_BSIZE;
        s = (enc->ts[WIDTH] + MAX_BSIZE - 1) / MAX_BSIZE;

        for (i = enc->quadtree_depth - 1; i >= 0; i--) {
            enc->qtmap[i] = (char ****) alloc_4d_array(v, u, t, s, sizeof(char));

            v <<= 1;
            u <<= 1;
            t <<= 1;
            s <<= 1;
        }
    }

    //Initiation of the reference offset
    enc->intra_roff = init_ref_offset(enc->vu, enc->ts, INTRA_PRED, enc->prd_order);

    // Keeps the weights used for the residue encoding context
    enc->intra_ctx_weight = init_ctx_weight(INTRA_PRED, enc->prd_order, enc->delta);

    if (enc->mi_prd_order[UP] > 0) {
        enc->mi_up_roff = init_ref_offset(enc->vu, enc->ts, MI_UP_PRED, enc->mi_prd_order[UP]);
        enc->mi_up_ctx_weight = init_ctx_weight(MI_UP_PRED, enc->mi_prd_order[UP], enc->delta);

        enc->mi_borders_ctx_weight = init_ctx_weight(MI_BORDERS_PRED, enc->mi_prd_order[UP], enc->delta);
    }

    if (enc->mi_prd_order[LEFT] > 0) {
        enc->mi_left_roff = init_ref_offset(enc->vu, enc->ts, MI_LEFT_PRED, enc->mi_prd_order[LEFT]);
        enc->mi_left_ctx_weight = init_ctx_weight(MI_LEFT_PRED, enc->mi_prd_order[LEFT], enc->delta);

        if (enc->mi_prd_order[LEFT] > enc->mi_prd_order[UP]) {
            enc->mi_borders_ctx_weight = init_ctx_weight(MI_BORDERS_PRED, enc->mi_prd_order[LEFT], enc->delta);
        }
    }

    if (enc->mi_prd_order[LDIAG] > 0) {
        enc->mi_ldiag_roff = init_ref_offset(enc->vu, enc->ts, MI_LDIAG_PRED, enc->mi_prd_order[LDIAG]);
        enc->mi_ldiag_ctx_weight = init_ctx_weight(MI_LDIAG_PRED, enc->mi_prd_order[LDIAG], enc->delta);

        if (enc->mi_prd_order[LDIAG] > enc->mi_prd_order[UP] && enc->mi_prd_order[LDIAG] > enc->mi_prd_order[LEFT]) {
            enc->mi_borders_ctx_weight = init_ctx_weight(MI_BORDERS_PRED, enc->mi_prd_order[LDIAG], enc->delta);
        }
    }

    if (enc->mi_prd_order[RDIAG] > 0) {
        enc->mi_rdiag_roff = init_ref_offset(enc->vu, enc->ts, MI_RDIAG_PRED, enc->mi_prd_order[RDIAG]);
        enc->mi_rdiag_ctx_weight = init_ctx_weight(MI_RDIAG_PRED, enc->mi_prd_order[RDIAG], enc->delta);

        if (enc->mi_prd_order[RDIAG] > enc->mi_prd_order[UP] && enc->mi_prd_order[RDIAG] > enc->mi_prd_order[LEFT] &&
            enc->mi_prd_order[RDIAG] > enc->mi_prd_order[LDIAG]) {
            enc->mi_borders_ctx_weight = init_ctx_weight(MI_BORDERS_PRED, enc->mi_prd_order[RDIAG], enc->delta);
        }
    }

    // Class and group arrays
    enc->class = (char ****) alloc_4d_array(enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH], sizeof(char));
    enc->group = (char ****) alloc_4d_array(enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH], sizeof(char));

    // Copy original images to encoder struct, initialize group array
    for (v = 0; v < enc->vu[HEIGHT]; v++) {
        for (u = 0; u < enc->vu[WIDTH]; u++) {
            for (t = 0; t < enc->ts[HEIGHT]; t++) {
                for (s = 0; s < enc->ts[WIDTH]; s++) {
                    enc->group[v][u][t][s] = 0;
                    enc->org[v][u][t][s] = lf->val[v][u][t][s];
                }
            }
        }
    }

    // Auxiliary values
    enc->org[enc->vu[HEIGHT]][0][0][0] = (enc->maxval + 1) >> 1;
    enc->err[enc->vu[HEIGHT]][0][0][0] = (enc->maxval + 1) >> 2;

    // Table used for the quantization of the context variable
    enc->uquant = (char **) alloc_2d_array(enc->num_class, MAX_UPARA + 1, sizeof(char));

    for (i = 0; i < enc->num_class; i++) {
        for (j = 0; j <= MAX_UPARA; j++) {
            enc->uquant[i][j] = (char) (enc->num_group - 1);
        }
    }

    // Table used for the conversion of the error
    if (enc->depth <= 8) {
        enc->econv = (int **) alloc_2d_array(enc->maxval + 1, (enc->maxval << 1) + 1, sizeof(int));
    }

    // Structure used to convert the prediction to a pointer which indicates the position in the probability vector structure of the prediction error
    enc->bconv = (img_t *) alloc_mem((enc->maxprd + 1) * sizeof(img_t));
    // Structure used to fine tune the probability value, given the probability model accuracy
    enc->fconv = (img_t *) alloc_mem((enc->maxprd + 1) * sizeof(img_t));
    // List of pointer for the probability model for each group
    enc->pmlist = (PMODEL **) alloc_mem(enc->num_group * sizeof(PMODEL *));
    // Probability model structure used for the side information (classes, thresholds, coefficients)
    enc->spm.freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
    enc->spm.cumfreq = &(enc->spm.freq[MAX_SYMBOL]);

    // Huffman coding
    if (enc->f_huffman == 1) {
        enc->sigma = sigma_h;
    }
    else {
        enc->sigma = sigma_a;
    }

    // Initialize probability models
    enc->pmodels = init_pmodels(enc->num_group, enc->num_pmodel, enc->pm_accuracy, NULL, enc->sigma, enc->maxval + 1);

    // Huffman coding
    if (enc->f_huffman == 1) {
        enc->vlcs = init_vlcs(enc->pmodels, enc->num_group, enc->num_pmodel);
    }

    // Buffer used for context determination for class encoding
    enc->mtfbuf = (int *) alloc_mem(enc->num_class * sizeof(int));
    // Structure that indicates the context used for arithmetic encoding of the coefficients of the prediction filters
    enc->coef_m = (int *) alloc_mem((enc->full_prd_order) * sizeof(int));

    for (i = 0; i < enc->full_prd_order; i++) {
        enc->coef_m[i] = 0;
    }

    // Structure used to keep the cost of the coefficients
    enc->coef_cost = (cost_t **) alloc_2d_array(16, MAX_COEF + 1, sizeof(cost_t));

    for (i = 0; i < 16; i++) {
#ifdef OPT_SIDEINFO
        if (enc->f_huffman == 1) {
            for (j = 0; j <= MAX_COEF; j++) {
                enc->coef_cost[i][j] = ((j >> i) + i + 1);
                if (j > 0) enc->coef_cost[i][j] += 1.0;
            }
        }
        else {
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

    // Structure used to keep the cost of the thresholds
    enc->th_cost = (cost_t *) alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
    for (i = 0; i < MAX_UPARA + 2; i++) {
        enc->th_cost[i] = 0;
    }

    // Array with the cost of each class
    enc->class_cost = (cost_t *) alloc_mem(enc->num_class * sizeof(cost_t));

    c = log((double) enc->num_class) / log(2.0);

    for (i = 0; i < enc->num_class; i++) {
        enc->class_cost[i] = c;
    }
    for (i = 0; i < (QUADTREE_DEPTH << 3); i++) {
        enc->qtflag_cost[i] = 1.0;
    }

    return (enc);
}

/*---------------------------- free_encoder ------------------------*
 |  Function free_encoder
 |
 |  Purpose:  Frees the memory used by the encoder structure
 |
 |  Parameters:
 |      enc				--> Encoder structure (IN)
 |
 |  Returns:  void
 *-------------------------------------------------------------------*/
void free_encoder(ENCODER *enc) {
    int x, i, j, k;
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
    if (enc->depth <= 8) {
        free(enc->econv);
    }
    free(enc->bconv);
    free(enc->fconv);
    free(enc->pmlist);
    free(enc->spm.freq);
    free(enc->mtfbuf);
    free(enc->coef_m);
    free(enc->coef_cost);
    free(enc->th_cost);
    free(enc->class_cost);

    free(enc->intra_ctx_weight);

    if(enc->mi_prd_order[UP] > 0) free(enc->mi_up_ctx_weight);
    if(enc->mi_prd_order[LEFT] > 0) free(enc->mi_left_ctx_weight);
    if(enc->mi_prd_order[LDIAG] > 0) free(enc->mi_ldiag_ctx_weight);
    if(enc->mi_prd_order[RDIAG] > 0) free(enc->mi_rdiag_ctx_weight);
    if(enc->mi_prd_order[UP] > 0 || enc->mi_prd_order[LEFT] > 0 || enc->mi_prd_order[LDIAG] > 0 || enc->mi_prd_order[RDIAG] > 0) free(enc->mi_borders_ctx_weight);

    if (enc->f_huffman == 0) {
        free(enc->rc);
    }

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

    free_ref_offset(enc->vu, enc->ts, INTRA_PRED, enc->prd_order, enc->intra_roff);

    if (enc->mi_prd_order[UP] > 0) free_ref_offset(enc->vu, enc->ts, MI_UP_PRED, enc->mi_prd_order[UP], enc->mi_up_roff);
    if (enc->mi_prd_order[LEFT] > 0) free_ref_offset(enc->vu, enc->ts, MI_LEFT_PRED, enc->mi_prd_order[LEFT], enc->mi_left_roff);
    if (enc->mi_prd_order[LDIAG] > 0) free_ref_offset(enc->vu, enc->ts, MI_LDIAG_PRED, enc->mi_prd_order[LDIAG], enc->mi_ldiag_roff);
    if (enc->mi_prd_order[RDIAG] > 0) free_ref_offset(enc->vu, enc->ts, MI_RDIAG_PRED, enc->mi_prd_order[RDIAG], enc->mi_rdiag_roff);

    free(enc->debug_path);

    free(enc);
}

/*------------------------------ init_class -------------------------*
 |  Function init_class
 |
 |  Purpose:  Sets the class number of each pixel in a given frame
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |
 |  Returns:  void
 *-------------------------------------------------------------------*/
void init_class(ENCODER *enc) {
    int k, v, u, t, s, g, h, i, j, z, cl, sum, num_block;
    int *var, *tmp, **ptr;

    int vu[2], ts[2], bsize_vu[2], bsize_ts[2];

    // Calculates the next multiple of BASE_BSIZE to correctly calculate the number of blocks
    vu[HEIGHT] = (int) ceil((double) enc->vu[HEIGHT] / BASE_BSIZE) * BASE_BSIZE;
    vu[WIDTH] = (int) ceil((double) enc->vu[WIDTH] / BASE_BSIZE) * BASE_BSIZE;

    ts[HEIGHT] = (int) ceil((double) enc->ts[HEIGHT] / BASE_BSIZE) * BASE_BSIZE;
    ts[WIDTH] = (int) ceil((double) enc->ts[WIDTH] / BASE_BSIZE) * BASE_BSIZE;

    // Number of blocks in the frame
    num_block = vu[HEIGHT] * vu[WIDTH] * ts[HEIGHT] * ts[WIDTH] / (BASE_BSIZE * BASE_BSIZE * BASE_BSIZE * BASE_BSIZE);

    var = (int *) alloc_mem(num_block * sizeof(int));
    ptr = (int **) alloc_mem(num_block * sizeof(int *));

    // Variance calculation for each block
    for (k = 0; k < num_block; k++) {
        // Gives the position of the top right pixel of each block
        v = (k / ((vu[WIDTH] * ts[HEIGHT] * ts[WIDTH]) / (BASE_BSIZE * BASE_BSIZE * BASE_BSIZE))) * BASE_BSIZE;
        u = ((k % ((vu[WIDTH] * ts[HEIGHT] * ts[WIDTH]) / (BASE_BSIZE * BASE_BSIZE * BASE_BSIZE))) / ((ts[HEIGHT] * ts[WIDTH]) / (BASE_BSIZE * BASE_BSIZE))) * BASE_BSIZE;

        t = ((k % ((ts[HEIGHT] * ts[WIDTH]) / (BASE_BSIZE * BASE_BSIZE))) / (ts[WIDTH] / BASE_BSIZE)) * BASE_BSIZE;
        s = (k % (ts[WIDTH] / BASE_BSIZE)) * BASE_BSIZE;

        var[k] = sum = 0;

        // Check the correct limits due to the image not being multiple of BASE_BSIZE
        bsize_vu[HEIGHT] = (v + BASE_BSIZE > enc->vu[HEIGHT]) ? enc->vu[HEIGHT] - v : BASE_BSIZE;
        bsize_vu[WIDTH] = (u + BASE_BSIZE > enc->vu[WIDTH]) ? enc->vu[WIDTH] - u : BASE_BSIZE;

        bsize_ts[HEIGHT] = (t + BASE_BSIZE > enc->ts[HEIGHT]) ? enc->ts[HEIGHT] - t : BASE_BSIZE;
        bsize_ts[WIDTH] = (s + BASE_BSIZE > enc->ts[WIDTH]) ? enc->ts[WIDTH] - s : BASE_BSIZE;

        // Run each pixel in a block
        for (g = 0; g < bsize_vu[HEIGHT]; g++) {
            for (h = 0; h < bsize_vu[WIDTH]; h++) {
                for (i = 0; i < bsize_ts[HEIGHT]; i++) {
                    for (j = 0; j < bsize_ts[WIDTH]; j++) {
                        z = enc->org[v + g][u + h][t + i][s + j];
                        sum += z;
                        var[k] += z * z;
                    }
                }
            }
        }

        // Final result of the variance for one block
        var[k] -= sum * sum / (bsize_vu[HEIGHT] * bsize_vu[WIDTH] * bsize_vu[HEIGHT] * bsize_ts[WIDTH]);
        ptr[k] = &(var[k]);
    }

    // Variance sorting, from lowest to highest
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
        // Calculates the number of the class for a block, given its position (1, 2, ...)
        cl = (k * enc->num_class) / num_block;
        // Determines the new position of a block after being sorted
        z = (int) (ptr[k] - var);

        // Gives the position of the top right pixel of each of the sorted blocks
        v = (z / ((vu[WIDTH] * ts[HEIGHT] * ts[WIDTH]) / (BASE_BSIZE * BASE_BSIZE * BASE_BSIZE))) * BASE_BSIZE;
        u = ((z % ((vu[WIDTH] * ts[HEIGHT] * ts[WIDTH]) / (BASE_BSIZE * BASE_BSIZE * BASE_BSIZE))) / ((ts[HEIGHT] * ts[WIDTH]) / (BASE_BSIZE * BASE_BSIZE))) * BASE_BSIZE;

        t = ((z % ((ts[HEIGHT] * ts[WIDTH]) / (BASE_BSIZE * BASE_BSIZE))) / (ts[WIDTH] / BASE_BSIZE)) * BASE_BSIZE;
        s = (z % (ts[WIDTH] / BASE_BSIZE)) * BASE_BSIZE;

        // Check the correct limits due to the image not being multiple of BASE_BSIZE
	bsize_vu[HEIGHT] = (v + BASE_BSIZE > enc->vu[HEIGHT]) ? enc->vu[HEIGHT] - v : BASE_BSIZE;
	bsize_vu[WIDTH] = (u + BASE_BSIZE > enc->vu[WIDTH]) ? enc->vu[WIDTH] - u : BASE_BSIZE;

	bsize_ts[HEIGHT] = (t + BASE_BSIZE > enc->ts[HEIGHT]) ? enc->ts[HEIGHT] - t : BASE_BSIZE;
	bsize_ts[WIDTH] = (s + BASE_BSIZE > enc->ts[WIDTH]) ? enc->ts[WIDTH] - s : BASE_BSIZE;

        // Sets the class number for each pixel
        for (g = 0; g < bsize_vu[HEIGHT]; g++) {
            for (h = 0; h < bsize_vu[WIDTH]; h++) {
                for (i = 0; i < bsize_ts[HEIGHT]; i++) {
                    for (j = 0; j < bsize_ts[WIDTH]; j++) {
                        enc->class[v + g][u + h][t + i][s + j] = (char) cl;
                    }
                }
            }
        }
    }

    free(ptr);
    free(var);
}

// Function to fill econv table
int error_conversion(int org, int prd, int maxval, int mode) {
    if (mode == 0) {
        int aux = (org << 1) - prd;

        return (aux > 0) ? (aux - 1) : (-aux);
    }
    else {
        int aux = (prd + 1) >> 1;
        return e2E(org - aux, aux, prd & 1, maxval);
    }
}

/*------------------------------ set_cost_model -------------------------*
 |  Function set_cost_model
 |
 |  Purpose:  Sets the cost model
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |		f_mmse			--> Sets the use of MMSE (IN)
 |
 |  Returns:  void
 *-------------------------------------------------------------------*/
void set_cost_model(ENCODER *enc, int f_mmse) {
    int gr, i, j, k;
    double a, b, c, var;
    PMODEL *pm;

    if (enc->depth <= 8) {
        for (i = 0; i <= enc->maxval; i++) {
            for (j = 0; j <= (enc->maxval << 1); j++) {
                enc->econv[i][j] = error_conversion(i, j, enc->maxval, enc->etype);
            }
        }
    }

    enc->encval = enc->err;

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
            c = (double) k * 0.5 + 0.25;
            pm->cost[k] = (float) (a + b * c * c);
        }

        pm->subcost[0] = 0.0;
    }

    for (k = 0; k <= enc->maxprd; k++) {
        enc->bconv[k] = 0;
        enc->fconv[k] = 0;
    }
}

/*------------------------------ set_cost_rate -------------------------*
 |  Function set_cost_rate
 |
 |  Purpose:  Sets the cost rate
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |
 |  Returns:  void
 *-------------------------------------------------------------------*/
void set_cost_rate(ENCODER *enc) {
    int gr, k, i, j, mask, shift, num_spm;
    double a, c;
    PMODEL *pm;

    if (enc->pm_accuracy < 0) {
        enc->etype = 1;

        if (enc->depth <= 8) {
            for (i = 0; i <= enc->maxval; i++) {
                for (j = 0; j <= (enc->maxval << 1); j++) {
                    enc->econv[i][j] = error_conversion(i, j, enc->maxval, enc->etype);
                }
            }
        }
    }

    if (enc->pm_accuracy < 0) {
        num_spm = 1;
    }
    else {
        enc->encval = enc->org;
        mask = (1 << enc->pm_accuracy) - 1;
        shift = enc->coef_precision - enc->pm_accuracy;

        for (k = 0; k <= enc->maxprd; k++) {
            i = (enc->maxprd - k + (1 << shift) / 2) >> shift;
            enc->fconv[k] = (img_t) (i & mask);
            enc->bconv[k] = (img_t) (i >> enc->pm_accuracy);
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
                    pm->cost[k] = (float) (-a * log(pm->freq[k]));
                }

                c = pm->cumfreq[enc->maxval + 1];
                pm->subcost[0] = (float) (a * log(c));
            }
            else {
                for (j = 0; j < num_spm; j++) {
                    for (k = 0; k < pm->size; k++) {
                        pm->cost[k] = (float) (-a * log(pm->freq[k]));
                    }

                    for (k = 0; k <= enc->maxval; k++) {
                        c = pm->cumfreq[k + enc->maxval + 1] - pm->cumfreq[k];
                        pm->subcost[k] = (float) (a * log(c));
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
 |  Purpose: Calculates the prediction of each pixel
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |		tlv				--> Top left V boundary (IN)
 |		tlu				--> Top left U boundary (IN)
 |		tlt				--> Top left T boundary (IN)
 |		tls				--> Top left S boundary (IN)
 |		brv				--> Bottom right V boundary (IN)
 |		bru				--> Bottom right U boundary (IN)
 |		brt				--> Bottom right T boundary (IN)
 |		brs				--> Bottom right S boundary (IN)
 |
 |  Returns:  void
 *-------------------------------------------------------------------*/
void predict_region(ENCODER *enc, int tlv, int tlu, int tlt, int tls, int brv, int bru, int brt, int brs) {
    int v, u, t, s, k, cl, prd, org;
    int *coef_p;
    int *prd_p;

    int *intra_roff_p = NULL, **intra_roff_pp = NULL;

    int *mi_up_roff_p = NULL, **mi_up_roff_pp = NULL;
    int *mi_left_roff_p = NULL, **mi_left_roff_pp = NULL;
    int *mi_ldiag_roff_p = NULL, **mi_ldiag_roff_pp = NULL;
    int *mi_rdiag_roff_p = NULL, **mi_rdiag_roff_pp = NULL;

    int *err_p, *org_p;
    char *class_p;

    //Runs all the pixels in a frame (due to the used boundaries)
    for (v = tlv; v < brv; v++) {
        for (u = tlu; u < bru; u++) {
            for (t = tlt; t < brt; t++) {
                class_p = &enc->class[v][u][t][tls];
                org_p = &enc->org[v][u][t][tls];
                err_p = &enc->err[v][u][t][tls];
                prd_p = &enc->prd[v][u][t][tls];

                intra_roff_pp = &enc->intra_roff[v][u][t][tls];

                if (enc->mi_prd_order[UP] > 0) mi_up_roff_pp = &enc->mi_up_roff[v][u][t][tls];
                if (enc->mi_prd_order[LEFT] > 0) mi_left_roff_pp = &enc->mi_left_roff[v][u][t][tls];
                if (enc->mi_prd_order[LDIAG] > 0) mi_ldiag_roff_pp = &enc->mi_ldiag_roff[v][u][t][tls];
                if (enc->mi_prd_order[RDIAG] > 0) mi_rdiag_roff_pp = &enc->mi_rdiag_roff[v][u][t][tls];

                for (s = tls; s < brs; s++) {
                    cl = *class_p++;
                    coef_p = enc->predictor[cl];

                    intra_roff_p = *intra_roff_pp++;

                    if (enc->mi_prd_order[UP] > 0) mi_up_roff_p = *mi_up_roff_pp++;
                    if (enc->mi_prd_order[LEFT] > 0) mi_left_roff_p = *mi_left_roff_pp++;
                    if (enc->mi_prd_order[LDIAG] > 0) mi_ldiag_roff_p = *mi_ldiag_roff_pp++;
                    if (enc->mi_prd_order[RDIAG] > 0) mi_rdiag_roff_p = *mi_rdiag_roff_pp++;

                    prd = 0;

                    for (k = 0; k < enc->prd_order; k++) {
                        prd += org_p[*intra_roff_p++] * (*coef_p++);
                    }

                    for (k = 0; k < enc->mi_prd_order[UP]; k++) {
                        prd += org_p[*mi_up_roff_p++] * (*coef_p++);
                    }

                    for (k = 0; k < enc->mi_prd_order[LEFT]; k++) {
                        prd += org_p[*mi_left_roff_p++] * (*coef_p++);
                    }

                    for (k = 0; k < enc->mi_prd_order[LDIAG]; k++) {
                        prd += org_p[*mi_ldiag_roff_p++] * (*coef_p++);
                    }

                    for (k = 0; k < enc->mi_prd_order[RDIAG]; k++) {
                        prd += org_p[*mi_rdiag_roff_p++] * (*coef_p++);
                    }

                    org = *org_p++;
                    *prd_p++ = prd;

                    //Boundaries of the predicted values
                    if (prd < 0) prd = 0;
                    else if (prd > enc->maxprd) prd = enc->maxprd;

                    prd >>= (enc->coef_precision - 1);

                    if (enc->depth <= 8) {
                        *err_p++ = enc->econv[org][prd];
                    }
                    else {
                        *err_p++ = error_conversion(org, prd, enc->maxval, enc->etype);
                    }
                }
            }
        }
    }
}

/*------------------------------ calc_uenc --------------------------*
 |  Function calc_uenc
 |
 |  Purpose: Calculates the context for the residue encoding,
 |			 without the quantization
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |		v				--> V position (IN)
 |		u				--> U position (IN)
 |		t				--> T position (IN)
 |		s				--> S position (IN)
 |
 |  Returns:  int		--> Returns the context value
 *-------------------------------------------------------------------*/
int calc_uenc(ENCODER *enc, int v, int u, int t, int s) {
    int uenc, k, *err_p;
    int *intra_wt_p;
    int *mi_up_wt_p = NULL, *mi_left_wt_p = NULL, *mi_ldiag_wt_p = NULL, *mi_rdiag_wt_p = NULL;
    int *intra_roff_p = NULL;
    int *mi_up_roff_p = NULL, *mi_left_roff_p = NULL, *mi_ldiag_roff_p = NULL, *mi_rdiag_roff_p = NULL;

    intra_roff_p = enc->intra_roff[v][u][t][s];
    intra_wt_p = enc->intra_ctx_weight;

    if (enc->mi_prd_order[UP] > 0) {
        mi_up_roff_p = enc->mi_up_roff[v][u][t][s];

        mi_up_wt_p = (v < 1) ? enc->mi_borders_ctx_weight : enc->mi_up_ctx_weight;
    }

    if (enc->mi_prd_order[LEFT] > 0) {
        mi_left_roff_p = enc->mi_left_roff[v][u][t][s];

        mi_left_wt_p = (u < 1) ? enc->mi_borders_ctx_weight : enc->mi_left_ctx_weight;
    }

    if (enc->mi_prd_order[LDIAG] > 0) {
        mi_ldiag_roff_p = enc->mi_ldiag_roff[v][u][t][s];

        mi_ldiag_wt_p = (v < 1 || u < 1) ? enc->mi_borders_ctx_weight : enc->mi_ldiag_ctx_weight;
    }

    if (enc->mi_prd_order[RDIAG] > 0) {
        mi_rdiag_roff_p = enc->mi_rdiag_roff[v][u][t][s];
        mi_rdiag_wt_p = (v < 1 || u > enc->vu[WIDTH] - 1) ? enc->mi_borders_ctx_weight : enc->mi_rdiag_ctx_weight;
    }

    err_p = &enc->err[v][u][t][s];

    uenc = 0;

    for (k = 0; k < enc->prd_order; k++) {
        uenc += err_p[*intra_roff_p++] * (*intra_wt_p++);
    }

    for (k = 0; k < enc->mi_prd_order[UP]; k++) {
        uenc += err_p[*mi_up_roff_p++] * (*mi_up_wt_p++);
    }

    for (k = 0; k < enc->mi_prd_order[LEFT]; k++) {
        uenc += err_p[*mi_left_roff_p++] * (*mi_left_wt_p++);
    }

    for (k = 0; k < enc->mi_prd_order[LDIAG]; k++) {
        uenc += err_p[*mi_ldiag_roff_p++] * (*mi_ldiag_wt_p++);
    }

    for (k = 0; k < enc->mi_prd_order[RDIAG]; k++) {
        uenc += err_p[*mi_rdiag_roff_p++] * (*mi_rdiag_wt_p++);
    }

    uenc >>= 6;

    if (uenc > MAX_UPARA) uenc = MAX_UPARA;

    return (uenc);
}

/*------------------------------ calc_cost --------------------------*
 |  Function calc_cost
 |
 |  Purpose: Iterates through all the pixel positions and sums up
 |			 the cost for encoding each prediction residue
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |		tlv				--> Top left V boundary (IN)
 |		tlu				--> Top left U boundary (IN)
 |		tlt				--> Top left T boundary (IN)
 |		tls				--> Top left S boundary (IN)
 |		brv				--> Bottom right V boundary (IN)
 |		bru				--> Bottom right U boundary (IN)
 |		brt				--> Bottom right T boundary (IN)
 |		brs				--> Bottom right S boundary (IN)
 |
 |  Returns:  cost_t	--> Returns the cost
 *-------------------------------------------------------------------*/
cost_t calc_cost(ENCODER *enc, int tlv, int tlu, int tlt, int tls, int brv, int bru, int brt, int brs) {
    cost_t cost;
    int v, u, t, s, uenc, cl, gr, prd, e, base, frac;
    int *upara_p, *prd_p, *encval_p;
    char *class_p, *group_p;
    PMODEL *pm;

    if (brv > enc->vu[HEIGHT]) brv = enc->vu[HEIGHT];
    if (bru > enc->vu[WIDTH]) bru = enc->vu[WIDTH];
    if (brt > enc->ts[HEIGHT]) brt = enc->ts[HEIGHT];
    if (tls < 0) tls = 0;
    if (brs > enc->ts[WIDTH]) brs = enc->ts[WIDTH];

    cost = 0;

    for (v = tlv; v < brv; v++) {
        for (u = tlu; u < bru; u++) {
            for (t = tlt; t < brt; t++) {
                class_p = &enc->class[v][u][t][tls];
                group_p = &enc->group[v][u][t][tls];
                upara_p = &enc->upara[v][u][t][tls];
                encval_p = &enc->encval[v][u][t][tls];
                prd_p = &enc->prd[v][u][t][tls];

                for (s = tls; s < brs; s++) {
                    cl = *class_p++;
                    *upara_p++ = uenc = calc_uenc(enc, v, u, t, s);
                    gr = enc->uquant[cl][uenc];
                    *group_p++ = (char) gr;
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
        }
    }

    return (cost);
}

/*-------------------------- design_predictor -----------------------*
 |  Function design_predictor
 |
 |  Purpose:  Designs the predictor coefficients for each class
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |		f_mmse			--> Coding type (Huffman/Arithmetic) (IN)
 |
 |  Returns:  cost_t --> cost of the designed predictor
 *-------------------------------------------------------------------*/
cost_t design_predictor(ENCODER *enc, int f_mmse) {
    double **mat, *weight, w, e, d, pivot;
    int v, u, t, s, i, j, k, cl, gr, pivpos, *index, *org_p;
    int *roff_p = NULL;

    mat = (double **) alloc_2d_array(enc->full_prd_order, enc->full_prd_order + 1, sizeof(double));
    index = (int *) alloc_mem(sizeof(int) * enc->full_prd_order);
    weight = (double *) alloc_mem(sizeof(double) * enc->num_group);

    // Weight choice
    for (gr = 0; gr < enc->num_group; gr++) {
        if (f_mmse) {
            weight[gr] = 1.0;
        }
        else {
            weight[gr] = 1.0 / (enc->sigma[gr] * enc->sigma[gr]);
        }
    }

    // Cycle that runs each class in a given frame
    for (cl = 0; cl < enc->num_class; cl++) {
        // Variable mat initialization
        for (i = 0; i < enc->full_prd_order; i++) {
            for (j = 0; j <= enc->full_prd_order; j++) {
                mat[i][j] = 0.0;
            }
        }

        // Run each pixel
        for (v = 0; v < enc->vu[HEIGHT]; v++) {
            for (u = 0; u < enc->vu[WIDTH]; u++) {
                for (t = 0; t < enc->ts[HEIGHT]; t++) {
                    for (s = 0; s < enc->ts[WIDTH]; s++) {
                        if (enc->class[v][u][t][s] != cl) { //Check if the pixel is member of the current class
                            s += BASE_BSIZE - 1;
                            continue;
                        }

                        gr = enc->group[v][u][t][s];
                        org_p = &enc->org[v][u][t][s];

                        roff_p = (int *) alloc_mem(sizeof(int) * enc->full_prd_order);

                        for (i = 0; i < enc->prd_order; i++) {
                            roff_p[i] = enc->intra_roff[v][u][t][s][i];
                        }

                        for (i = 0; i < enc->mi_prd_order[UP]; i++) {
                            roff_p[i + enc->prd_order] = enc->mi_up_roff[v][u][t][s][i];
                        }

                        for (i = 0; i < enc->mi_prd_order[LEFT]; i++) {
                            roff_p[i + enc->prd_order + enc->mi_prd_order[UP]] = enc->mi_left_roff[v][u][t][s][i];
                        }

                        for (i = 0; i < enc->mi_prd_order[LDIAG]; i++) {
                            roff_p[i + enc->prd_order + enc->mi_prd_order[UP] +
                                   enc->mi_prd_order[LEFT]] = enc->mi_ldiag_roff[v][u][t][s][i];
                        }

                        for (i = 0; i < enc->mi_prd_order[RDIAG]; i++) {
                            roff_p[i + enc->prd_order + enc->mi_prd_order[UP] + enc->mi_prd_order[LEFT] +
                                   enc->mi_prd_order[LDIAG]] = enc->mi_rdiag_roff[v][u][t][s][i];
                        }

                        // Fills the matrix mat with the reference pixels values
                        for (i = 0; i < enc->full_prd_order; i++) {
                            w = weight[gr] * org_p[roff_p[i]];

                            for (j = i; j < enc->full_prd_order; j++) {
                                mat[i][j] += w * org_p[roff_p[j]];
                            }

                            // Results column with the value of the pixel to be predicted
                            mat[i][enc->full_prd_order] += w * org_p[0];
                        }

                        free(roff_p);
                    }
                }
            }
        }

        // Makes the bottom diagonal equal to the top diagonal
        for (i = 0; i < enc->full_prd_order; i++) {
            index[i] = i;

            for (j = 0; j < i; j++) {
                mat[i][j] = mat[j][i];
            }
        }

        // "Gaussian elimination"
        for (i = 0; i < enc->full_prd_order; i++) {
            pivpos = i;
            pivot = fabs(mat[index[i]][i]);

            // Sorts the matrix pivots
            for (k = i + 1; k < enc->full_prd_order; k++) {
                if (fabs(mat[index[k]][i]) > pivot) {
                    pivot = fabs(mat[index[k]][i]);
                    pivpos = k;
                }
            }

            // Change the indexes
            k = index[i];
            index[i] = index[pivpos];
            index[pivpos] = k;

            // If the pivot is not zero the actual Gaussian elimination is performed
            if (pivot > 1E-10) {
                d = mat[index[i]][i];

                for (j = i; j <= enc->full_prd_order; j++) {
                    mat[index[i]][j] /= d;
                }

                for (k = 0; k < enc->full_prd_order; k++) {
                    if (k == i) continue;

                    d = mat[index[k]][i];

                    for (j = i; j <= enc->full_prd_order; j++) {
                        mat[index[k]][j] -= d * mat[index[i]][j];
                    }
                }
            }
        }

        w = (1 << enc->coef_precision);
        e = 0.0;

        // Rounds the coefficients and stores the error
        for (i = 0; i < enc->full_prd_order; i++) {
            if (fabs(mat[index[i]][i]) > 1E-10) { //Checks if a line is not zero
                d = mat[index[i]][enc->full_prd_order] * w;
            }
            else {
                d = 0.0;
            }

            k = (int) d;
            if (k > d) k--;

            if (k < -MAX_COEF) {
                d = k = -MAX_COEF;
            }
            else if (k > MAX_COEF) {
                d = k = MAX_COEF;
            }

            enc->predictor[cl][i] = k; // Coefficient for a given reference in the predictor
            d -= k;
            e += d;
            mat[index[i]][enc->full_prd_order] = d; // Stores the error
        }

        /* minimize mean rounding errors */
        k = (int) (e + 0.5);
        for (; k > 0; k--) {
            d = 0;
            for (j = i = 0; i < enc->full_prd_order; i++) {
                if (mat[index[i]][enc->full_prd_order] > d) {
                    d = mat[index[i]][enc->full_prd_order];
                    j = i;
                }
            }

            if (enc->predictor[cl][j] < MAX_COEF) enc->predictor[cl][j]++;

            mat[index[j]][enc->full_prd_order] = 0;
        }
    }

    free(weight);
    free(index);
    free(mat);

    predict_region(enc, 0, 0, 0, 0, enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH]);

    return (calc_cost(enc, 0, 0, 0, 0, enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH]));
}

/*--------------------------- optimize_group ------------------------*
 |  Function optimize_group
 |
 |  Purpose: Optimize the thresholds that determine the group
 |			 for each residue
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN)
 |
 |  Returns:  cost_t	--> Returns the new cost of encoding
 *-------------------------------------------------------------------*/
cost_t optimize_group(ENCODER *enc) {
    cost_t cost, min_cost, **cbuf, *dpcost, *cbuf_p, *thc_p;
    int v, u, t, s, th1, th0, k, upara, cl, gr, prd, e, base, frac;
    int **trellis, *tre_p;
    PMODEL *pm, **pm_p;

    trellis = (int **) alloc_2d_array(enc->num_group, MAX_UPARA + 2, sizeof(int));
    dpcost = (cost_t *) alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
    cbuf = (cost_t **) alloc_2d_array(enc->num_group, MAX_UPARA + 2, sizeof(cost_t));
    thc_p = enc->th_cost;

    for (k = 0; k < MAX_UPARA + 2; k++) trellis[0][k] = 0;
    /* Dynamic programming */
    for (cl = 0; cl < enc->num_class; cl++) {
        for (gr = 0; gr < enc->num_group; gr++) {
            cbuf_p = cbuf[gr];

            for (upara = 0; upara < MAX_UPARA + 2; upara++) {
                cbuf_p[upara] = 0;
            }
        }

        for (v = 0; v < enc->vu[HEIGHT]; v++) {
            for (u = 0; u < enc->vu[WIDTH]; u++) {
                for (t = 0; t < enc->ts[HEIGHT]; t++) {
                    for (s = 0; s < enc->ts[WIDTH]; s++) {
                        if (enc->class[v][u][t][s] == cl) {
                            upara = enc->upara[v][u][t][s] + 1;
                            e = enc->encval[v][u][t][s];
                            prd = enc->prd[v][u][t][s];

                            if (prd < 0) prd = 0;
                            else if (prd > enc->maxprd) prd = enc->maxprd;

                            base = enc->bconv[prd];
                            frac = enc->fconv[prd];
                            pm_p = enc->pmlist;

                            for (gr = 0; gr < enc->num_group; gr++) {
                                pm = (*pm_p++) + frac;
                                cbuf[gr][upara] += pm->cost[base + e] + pm->subcost[base];
                            }
                        }
                    }
                }
            }
        }

        for (gr = 0; gr < enc->num_group; gr++) {
            cbuf_p = cbuf[gr];

            for (upara = 1; upara < MAX_UPARA + 2; upara++) {
                cbuf_p[upara] += cbuf_p[upara - 1];
            }
        }

        cbuf_p = cbuf[0];

        for (upara = 0; upara < MAX_UPARA + 2; upara++) {
            dpcost[upara] = cbuf_p[upara];
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
        upara = 0;

        for (gr = 0; gr < enc->num_group; gr++) {
            for (; upara < enc->th[cl][gr]; upara++) {
                enc->uquant[cl][upara] = (char) gr;
            }
        }
    }

    /* renew groups */
    cost = 0;
    pm_p = enc->pmlist;

    for (v = 0; v < enc->vu[HEIGHT]; v++) {
        for (u = 0; u < enc->vu[WIDTH]; u++) {
            for (t = 0; t < enc->ts[HEIGHT]; t++) {
                for (s = 0; s < enc->ts[WIDTH]; s++) {
                    cl = enc->class[v][u][t][s];
                    upara = enc->upara[v][u][t][s];
                    gr = enc->group[v][u][t][s] = enc->uquant[cl][upara];
                    e = enc->encval[v][u][t][s];
                    prd = enc->prd[v][u][t][s];

                    if (prd < 0) prd = 0;
                    else if (prd > enc->maxprd) prd = enc->maxprd;

                    base = enc->bconv[prd];
                    pm = pm_p[gr] + enc->fconv[prd];
                    cost += pm->cost[base + e] + pm->subcost[base];
                }
            }
        }
    }

    /* optimize probability models */
    if (enc->optimize_loop > 1 && enc->num_pmodel > 1) {
        if (enc->num_pmodel > MAX_UPARA + 2) {
            free(cbuf);
            cbuf = (cost_t **) alloc_2d_array(enc->num_group, enc->num_pmodel, sizeof(cost_t));
        }

        for (gr = 0; gr < enc->num_group; gr++) {
            for (k = 0; k < enc->num_pmodel; k++) {
                cbuf[gr][k] = 0;
            }
        }

        for (v = 0; v < enc->vu[HEIGHT]; v++) {
            for (u = 0; u < enc->vu[WIDTH]; u++) {
                for (t = 0; t < enc->ts[HEIGHT]; t++) {
                    for (s = 0; s < enc->ts[WIDTH]; s++) {
                        gr = enc->group[v][u][t][s];
                        e = enc->encval[v][u][t][s];
                        prd = enc->prd[v][u][t][s];

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
 |  Purpose: Calculates and stores the prediction and residue
 |			 for all classes for a given block
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN)
 |		prdbuf			--> Struct to store the prediction value (IN/OUT)
 |		errbuf          --> Struct to store the residue value (IN/OUT)
 |		tlv				--> Top left corner v position of the current block (IN)
 |		tlu             --> Top left corner u position of the current block (IN)
 |		tlt				--> Top left corner t position of the current block (IN)
 |		tls             --> Top left corner s position of the current block (IN)
 |		bufsize			--> Size of the block (IN)
 |
 |  Returns:  ENCODER* --> returns a encoder type structure
 *-------------------------------------------------------------------*/
void set_prdbuf(ENCODER *enc, int **prdbuf, int **errbuf, int tlv, int tlu, int tlt, int tls, int bufsize) {
    int v, u, t, s, brv, bru, brt, brs, cl, k, prd, *prdbuf_p, *errbuf_p, *coef_p;
    int buf_ptr, org, *org_p;
    int *intra_roff_p = NULL;
    int *mi_up_roff_p = NULL, *mi_left_roff_p = NULL, *mi_ldiag_roff_p = NULL, *mi_rdiag_roff_p = NULL;

    //Check the boundaries of the block (mainly for pictures not multiple of bufsize)
    brv = (tlv + bufsize < enc->vu[HEIGHT]) ? (tlv + bufsize) : enc->vu[HEIGHT];
    bru = (tlu + bufsize < enc->vu[WIDTH]) ? (tlu + bufsize) : enc->vu[WIDTH];
    brt = (tlt + bufsize < enc->ts[HEIGHT]) ? (tlt + bufsize) : enc->ts[HEIGHT];
    brs = (tls + bufsize < enc->ts[WIDTH]) ? (tls + bufsize) : enc->ts[WIDTH];

    // Run all classes
    for (cl = 0; cl < enc->num_class; cl++) {
        //Run all pixels in the given block
        for (v = tlv; v < brv; v++) {
            for (u = tlu; u < bru; u++) {
                buf_ptr = bufsize * bufsize * bufsize * (v % bufsize) + bufsize * bufsize * (u % bufsize) +
                          bufsize * (tlt % bufsize) + tls % bufsize;

                for (t = tlt; t < brt; t++) {
                    prdbuf_p = &prdbuf[cl][buf_ptr];
                    errbuf_p = &errbuf[cl][buf_ptr];
                    buf_ptr += bufsize;
                    org_p = &enc->org[v][u][t][tls];

                    for (s = tls; s < brs; s++) {
                        //If the pixel is of the current class just copy
                        if (cl == enc->class[v][u][t][s]) {
                            *prdbuf_p++ = enc->prd[v][u][t][s];
                            *errbuf_p++ = enc->err[v][u][t][s];
                            org_p++;
                        }
                        else { //If not do the calculations
                            intra_roff_p = enc->intra_roff[v][u][t][s];

                            if (enc->mi_prd_order[UP] > 0) mi_up_roff_p = enc->mi_up_roff[v][u][t][s];
                            if (enc->mi_prd_order[LEFT] > 0) mi_left_roff_p = enc->mi_left_roff[v][u][t][s];
                            if (enc->mi_prd_order[LDIAG] > 0) mi_ldiag_roff_p = enc->mi_ldiag_roff[v][u][t][s];
                            if (enc->mi_prd_order[RDIAG] > 0) mi_rdiag_roff_p = enc->mi_rdiag_roff[v][u][t][s];

                            coef_p = enc->predictor[cl];
                            prd = 0;

                            for (k = 0; k < enc->prd_order; k++) {
                                prd += org_p[*intra_roff_p++] * (*coef_p++);
                            }

                            for (k = 0; k < enc->mi_prd_order[UP]; k++) {
                                prd += org_p[*mi_up_roff_p++] * (*coef_p++);
                            }

                            for (k = 0; k < enc->mi_prd_order[LEFT]; k++) {
                                prd += org_p[*mi_left_roff_p++] * (*coef_p++);
                            }

                            for (k = 0; k < enc->mi_prd_order[LDIAG]; k++) {
                                prd += org_p[*mi_ldiag_roff_p++] * (*coef_p++);
                            }

                            for (k = 0; k < enc->mi_prd_order[RDIAG]; k++) {
                                prd += org_p[*mi_rdiag_roff_p++] * (*coef_p++);
                            }

                            org = *org_p++;
                            *prdbuf_p++ = prd;

                            if (prd < 0) prd = 0;
                            else if (prd > enc->maxprd) prd = enc->maxprd;

                            prd >>= (enc->coef_precision - 1);

                            if (enc->depth <= 8) {
                                *errbuf_p++ = enc->econv[org][prd];
                            }
                            else {
                                *errbuf_p++ = error_conversion(org, prd, enc->maxval, enc->etype);
                            }
                        }
                    }
                }
            }
        }
    }
}

/*---------------------------- find_class ---------------------------*
 |  Function find_class
 |
 |  Purpose: Tests all possible prediction modes for the block
 |			 and returns the class that resulted in the minimum cost
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN)
 |		prdbuf			--> Struct to store the prediction value (IN/OUT)
 |		errbuf          --> Struct to store the residue value (IN/OUT)
 |		tly				--> Top left corner y position of the current block (IN)
 |		tlx             --> Top left corner x position of the current block (IN)
 |		bry				--> Bottom right Y position of the current block (IN)
 |		brx				--> Bottom right X position of the current block (IN)
 |		bufsize			--> Size of the block (IN)
 |
 |  Returns:  int 		--> Returns the best class
 *-------------------------------------------------------------------*/
int find_class(ENCODER *enc, int **prdbuf, int **errbuf, int tlv, int tlu, int tlt, int tls, int brv, int bru, int brt,
               int brs, int bufsize) {
    cost_t cost, min_cost;
    int v, u, t, s, bufptr, cl, min_cl;
    char *class_p;
    int *prd_p, *prdbuf_p, *err_p, *errbuf_p;

    min_cost = 1E8;
    min_cl = 0;

    for (cl = 0; cl < enc->num_class; cl++) {
        cost = enc->class_cost[enc->mtfbuf[cl]];

        for (v = tlv; v < brv; v++) {
            for (u = tlu; u < bru; u++) {
                bufptr = bufsize * bufsize * bufsize * (v % bufsize) + bufsize * bufsize * (u % bufsize) +
                         bufsize * (tlt % bufsize) + tls % bufsize;

                for (t = tlt; t < brt; t++) {
                    class_p = &enc->class[v][u][t][tls];
                    prd_p = &enc->prd[v][u][t][tls];
                    prdbuf_p = &prdbuf[cl][bufptr];
                    err_p = &enc->err[v][u][t][tls];
                    errbuf_p = &errbuf[cl][bufptr];
                    bufptr += bufsize;

                    for (s = tls; s < brs; s++) {
                        *class_p++ = (char) cl;
                        *prd_p++ = *prdbuf_p++;
                        *err_p++ = *errbuf_p++;
                    }
                }
            }
        }

        cost += calc_cost(enc, tlv, tlu, tlt, tls, brv, bru, brt, brs);

        if (cost < min_cost) {
            min_cost = cost;
            min_cl = cl;
        }
    }

    for (v = tlv; v < brv; v++) {
        for (u = tlu; u < bru; u++) {
            bufptr = bufsize * bufsize * bufsize * (v % bufsize) + bufsize * bufsize * (u % bufsize) +
                     bufsize * (tlt % bufsize) + tls % bufsize;

            for (t = tlt; t < brt; t++) {
                class_p = &enc->class[v][u][t][tls];
                prd_p = &enc->prd[v][u][t][tls];
                prdbuf_p = &prdbuf[min_cl][bufptr];
                err_p = &enc->err[v][u][t][tls];
                errbuf_p = &errbuf[min_cl][bufptr];
                bufptr += bufsize;

                for (s = tls; s < brs; s++) {
                    *class_p++ = (char) min_cl;
                    *prd_p++ = *prdbuf_p++;
                    *err_p++ = *errbuf_p++;
                }
            }
        }
    }

    return (min_cl);
}

/*------------------------------ vbs_class --------------------------*
 |  Function vbs_class
 |
 |  Purpose: Performs the search for the best class based
 |			 on the residue cost
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN)
 |		prdbuf			--> Struct to store the prediction value (IN/OUT)
 |		errbuf          --> Struct to store the residue value (IN/OUT)
 |		tlv				--> Top left corner v position of the current block (IN)
 |		tlu             --> Top left corner u position of the current block (IN)
 |		tlt				--> Top left corner t position of the current block (IN)
 |		tls             --> Top left corner s position of the current block (IN)
 |		blksize			--> Size of the block (IN)
 |		width			--> Image width (IN)
 |		level			--> Level of the quadtree (IN)
 |
 |  Returns:  cost_t --> The class and partition cost
 *-------------------------------------------------------------------*/
cost_t vbs_class(ENCODER *enc, int **prdbuf, int **errbuf, int tlv, int tlu, int tlt, int tls, int blksize, int width,
                 int level) {
    int v, u, t, s, k, brv, bru, brt, brs, cl, bufsize, bufptr, ctx;
    int mtf_save[MAX_CLASS];
    char ****qtmap;
    cost_t cost1, cost2, qtcost;
    char *class_p;
    int *prd_p, *prdbuf_p, *err_p, *errbuf_p;

    // Sets the max block size depending on the current loop and if the quadtree is active
    if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
        bufsize = MAX_BSIZE;
    }
    else {
        bufsize = BASE_BSIZE;
    }

    // Check the boundaries of the block (mainly for pictures not multiple of bufsize)
    brv = (tlv + blksize < enc->vu[HEIGHT]) ? (tlv + blksize) : enc->vu[HEIGHT];
    bru = (tlu + blksize < enc->vu[WIDTH]) ? (tlu + blksize) : enc->vu[WIDTH];
    brt = (tlt + blksize < enc->ts[HEIGHT]) ? (tlt + blksize) : enc->ts[HEIGHT];
    brs = (tls + blksize < enc->ts[WIDTH]) ? (tls + blksize) : enc->ts[WIDTH];

    // Safeguard condition
    if (tlv >= brv || tlu >= bru || tlt >= brt || tls >= brs) return (0);

    // Auxiliary variable that stores class order
    for (k = 0; k < enc->num_class; k++) {
        mtf_save[k] = enc->mtfbuf[k];
    }

    // In the case of the first loop or if we are not using quadtree
    // Reorders classes to favor the modes chosen by the neighboring blocks.
    mtf_classlabel(enc->class, enc->mtfbuf, tlv, tlu, tlt, tls, blksize, width, enc->num_class);
    // Finds the class that returns the minimum cost for the block
    cl = find_class(enc, prdbuf, errbuf, tlv, tlu, tlt, tls, brv, bru, brt, brs, bufsize);
    // Return the class cost
    qtcost = enc->class_cost[enc->mtfbuf[cl]];

    // In the case of the second loop and if we are using quadtree
    // Check the best block partitioning
    if (level > 0) {
        /* context for quad-tree flag */
        ctx = 0;

        qtmap = enc->qtmap[level - 1];

        v = (tlv / MIN_BSIZE) >> level;
        u = (tlu / MIN_BSIZE) >> level;
        t = (tlt / MIN_BSIZE) >> level;
        s = (tls / MIN_BSIZE) >> level;

        if (v > 0 && qtmap[v - 1][u][t][s] == 1) ctx++;

        if (u > 0 && qtmap[v][u - 1][t][s] == 1) ctx++;

        if (t > 0) {
            if (qtmap[v][u][t - 1][s] == 1) ctx++;
            if (brs < width && qtmap[v][u][t - 1][s + 1] == 1) ctx++;
        }

        if (s > 0 && qtmap[v][u][t][s - 1] == 1) ctx++;

        ctx = ((level - 1) * 4 + ctx) << 1;

        /* Quad-tree partitioning */
        cost1 = calc_cost(enc, tlv, tlu, tlt, tls, brv, bru, brt, brs) + enc->class_cost[enc->mtfbuf[cl]] + enc->qtflag_cost[ctx];
        blksize >>= 1;

        for (k = 0; k < enc->num_class; k++) {
            enc->mtfbuf[k] = mtf_save[k];
        }

        qtcost = enc->qtflag_cost[ctx + 1];

        // v and u
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv, tlu, tlt, tls, blksize, width, level - 1);
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv, tlu, tlt, tls + blksize, blksize, width, level - 1);
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv, tlu, tlt + blksize, tls, blksize, width, level - 1);
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv, tlu, tlt + blksize, tls + blksize, blksize, brs, level - 1);
        // v and u + blksize
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv, tlu + blksize, tlt, tls, blksize, width, level - 1);
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv, tlu + blksize, tlt, tls + blksize, blksize, width, level - 1);
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv, tlu + blksize, tlt + blksize, tls, blksize, width, level - 1);
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv, tlu + blksize, tlt + blksize, tls + blksize, blksize, brs, level - 1);
        // v + blksize and u
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv + blksize, tlu, tlt, tls, blksize, width, level - 1);
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv + blksize, tlu, tlt, tls + blksize, blksize, width, level - 1);
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv + blksize, tlu, tlt + blksize, tls, blksize, width, level - 1);
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv + blksize, tlu, tlt + blksize, tls + blksize, blksize, brs, level - 1);
        // v + blksize and u + blksize
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv + blksize, tlu + blksize, tlt, tls, blksize, width, level - 1);
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv + blksize, tlu + blksize, tlt, tls + blksize, blksize, width, level - 1);
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv + blksize, tlu + blksize, tlt + blksize, tls, blksize, width, level - 1);
        qtcost += vbs_class(enc, prdbuf, errbuf, tlv + blksize, tlu + blksize, tlt + blksize, tls + blksize, blksize, brs, level - 1);

        cost2 = calc_cost(enc, tlv, tlu, tlt, tls, brv, bru, brt, brs) + qtcost;

        if (cost1 < cost2) {
            blksize <<= 1;

            for (k = 0; k < enc->num_class; k++) {
                enc->mtfbuf[k] = mtf_save[k];
            }

            mtf_classlabel(enc->class, enc->mtfbuf, tlv, tlu, tlt, tls, blksize, width, enc->num_class);
            qtcost = enc->class_cost[enc->mtfbuf[cl]] + enc->qtflag_cost[ctx];

            for (v = tlv; v < brv; v++) {
                for (u = tlu; u < bru; u++) {
                    bufptr = bufsize * bufsize * bufsize * (v % bufsize) + bufsize * bufsize * (u % bufsize) +
                             bufsize * (tlt % bufsize) + tls % bufsize;

                    for (t = tlt; t < brt; t++) {
                        class_p = &enc->class[v][u][t][tls];
                        prd_p = &enc->prd[v][u][t][tls];
                        prdbuf_p = &prdbuf[cl][bufptr];
                        err_p = &enc->err[v][u][t][tls];
                        errbuf_p = &errbuf[cl][bufptr];
                        bufptr += bufsize;

                        for (s = tls; s < brs; s++) {
                            *class_p++ = (char) cl;
                            *prd_p++ = *prdbuf_p++;
                            *err_p++ = *errbuf_p++;
                        }
                    }
                }
            }

            tlv = (tlv / MIN_BSIZE) >> level;
            tlu = (tlu / MIN_BSIZE) >> level;
            tlt = (tlt / MIN_BSIZE) >> level;
            tls = (tls / MIN_BSIZE) >> level;

            brv = tlv + 1;
            bru = tlu + 1;
            brt = tlt + 1;
            brs = tls + 1;

            for (; level > 0; level--) {
                qtmap = enc->qtmap[level - 1];

                for (v = tlv; v < brv; v++) {
                    for (u = tlu; u < bru; u++) {
                        for (t = tlt; t < brt; t++) {
                            for (s = tls; s < brs; s++) {
                                qtmap[v][u][t][s] = 0;
                            }
                        }
                    }
                }

                tlv <<= 1;
                tlu <<= 1;
                tlt <<= 1;
                tls <<= 1;

                brv <<= 1;
                bru <<= 1;
                brt <<= 1;
                brs <<= 1;
            }
        }
        else {
            qtmap[v][u][t][s] = 1;
        }
    }

    return (qtcost);
}

/*--------------------------- optimize_class ------------------------*
 |  Function optimize_class
 |
 |  Purpose: Optimizes the prediction class for each block
 |			 in the picture
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN/OUT)
 |
 |  Returns:  ENCODER* --> returns a encoder type structure
 *-------------------------------------------------------------------*/
cost_t optimize_class(ENCODER *enc) {
    int v, u, t, s, i, blksize, level;
    int **prdbuf, **errbuf;

    // Sets the max block size depending on the current loop
    // and if the quadtree is active
    if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
        level = enc->quadtree_depth;
        blksize = MAX_BSIZE; // 32
    }
    else {
        level = 0;
        blksize = BASE_BSIZE; // 8
    }

    // Initialize buffer with the class numbers
    for (i = 0; i < enc->num_class; i++) {
        enc->mtfbuf[i] = i;
    }

    // Auxiliary structures, prediction and residue buffer
    prdbuf = (int **) alloc_2d_array(enc->num_class, blksize * blksize * blksize * blksize, sizeof(int));
    errbuf = (int **) alloc_2d_array(enc->num_class, blksize * blksize * blksize * blksize, sizeof(int));

    // Cycle to run all picture blocks
    for (v = 0; v < enc->vu[HEIGHT]; v += blksize) {
        for (u = 0; u < enc->vu[WIDTH]; u += blksize) {
            for (t = 0; t < enc->ts[HEIGHT]; t += blksize) {
                for (s = 0; s < enc->ts[HEIGHT]; s += blksize) {
                    // Calculates and stores the prediction and residue for all classes for a given block
                    set_prdbuf(enc, prdbuf, errbuf, v, u, t, s, blksize);
                    // Determines the variable size block partition
                    vbs_class(enc, prdbuf, errbuf, v, u, t, s, blksize, enc->ts[WIDTH], level);
                }
            }
        }
    }

    // Free auxiliary pointers
    free(errbuf);
    free(prdbuf);

    // Returns cost
    return (calc_cost(enc, 0, 0, 0, 0, enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH]));
}

/*--------------------------- optimize_coef -------------------------*
 |  Function optimize_coef
 |
 |  Purpose: Optimizes the coefficients
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN/OUT)
 |		cl				--> Class number (IN)
 |		pos1			--> Coefficient position (IN)
 |		pos2			--> Coefficient position (IN)
 |
 |  Returns:  void
 *-------------------------------------------------------------------*/
void optimize_coef(ENCODER *enc, int cl, int pos1, int pos2) {
//#define SEARCH_RANGE 11
//#define SUBSEARCH_RANGE 3
//	cost_t cbuf[SEARCH_RANGE * SUBSEARCH_RANGE], *cbuf_p;
//	float *pmcost_p;
//	int i, j, k, x, y, df1, df2, base;
//	int prd, prd_f, shift, maxprd, *coef_p, *econv_p = NULL, *org_p;
//	int ***roff1_pp, ***roff2_pp, *roff1_p, *roff2_p;
//    char *class_p;
//	img_t *bconv_p, *fconv_p;
//	PMODEL *pm, *pm_p;
//
//	cbuf_p = cbuf;
//	coef_p = enc->predictor[cl];
//	k = 0;
//
//	for (i = 0; i < SEARCH_RANGE; i++) {
//		y = coef_p[pos1] + i - (SEARCH_RANGE >> 1);
//
//		if (y < 0) y = -y;
//
//		if (y > MAX_COEF) y = MAX_COEF;
//
//		for (j = 0; j < SUBSEARCH_RANGE; j++) {
//			x = coef_p[pos2] - (i - (SEARCH_RANGE >> 1)) - (j - (SUBSEARCH_RANGE >> 1));
//
//			if (x < 0) x = -x;
//
//			if (x > MAX_COEF) x = MAX_COEF;
//
//			cbuf_p[k++] = enc->coef_cost[enc->coef_m[pos1]][y] + enc->coef_cost[enc->coef_m[pos2]][x];
//		}
//	}
//
//	bconv_p = enc->bconv;
//	fconv_p = enc->fconv;
//	maxprd = enc->maxprd;
//	shift = enc->coef_precision - 1;
//
//	// TODO: check case where for instance mi_prd_order == 0
//    roff1_pp = (pos1 < enc->prd_order) ? enc->intra_roff : (pos1 < enc->prd_order + enc->mi_prd_order[UP]) ? enc->mi_up_roff : (pos1 < enc->prd_order + enc->mi_prd_order[UP] + enc->back_prd_order) ? enc->back_roff : enc->for_roff;
//    roff2_pp = (pos2 < enc->prd_order) ? enc->intra_roff : (pos2 < enc->prd_order + enc->mi_prd_order[UP]) ? enc->mi_up_roff : (pos2 < enc->prd_order + enc->mi_prd_order[UP] + enc->back_prd_order) ? enc->back_roff : enc->for_roff;
//
//	for (y = 0; y < enc->height; y++) {
//		class_p = enc->class[y];
//
//		for (x = 0; x < enc->width; x++) {
//			if (cl != *class_p++) continue;
//
//			prd = enc->prd[y][x];
//			org_p = &enc->org[1][y][x];
//			pm_p = enc->pmlist[(int)enc->group[y][x]];
//
//            roff1_p = roff1_pp[y][x];
//            roff2_p = roff2_pp[y][x];
//
//            df1 = org_p[roff1_p[pos1]];
//			df2 = org_p[roff2_p[pos2]];
//
//			prd_f = prd - (df1 - df2) * (SEARCH_RANGE >> 1) + df2 * (SUBSEARCH_RANGE >> 1);
//			cbuf_p = cbuf;
//
//			if (enc->pm_accuracy < 0) {
//				if (enc->depth <= 8) {
//					econv_p = enc->econv[*org_p];
//				}
//
//				pmcost_p = pm_p->cost;
//
//				for (i = 0; i < SEARCH_RANGE; i++) {
//					for (j = 0; j < SUBSEARCH_RANGE; j++) {
//						prd = prd_f;
//
//						if (prd < 0) prd = 0;
//						else if (prd > maxprd) prd = maxprd;
//
//						if (enc->depth <= 8) {
//							(*cbuf_p++) += pmcost_p[econv_p[prd >> shift]];
//						}
//						else {
//							(*cbuf_p++) += pmcost_p[error_conversion(*org_p, prd >> shift, enc->maxval, enc->etype)];
//						}
//
//						prd_f -= df2;
//					}
//
//					prd_f += df1 + df2 * (SUBSEARCH_RANGE - 1);
//				}
//			}
//			else {
//				for (i = 0; i < SEARCH_RANGE; i++) {
//					for (j = 0; j < SUBSEARCH_RANGE; j++) {
//						prd = prd_f;
//
//						if (prd < 0) prd = 0;
//						else if (prd > maxprd) prd = maxprd;
//
//						base = bconv_p[prd];
//						pm = pm_p + fconv_p[prd];
//						(*cbuf_p++) += pm->cost[*org_p + base] + pm->subcost[base];
//						prd_f -= df2;
//					}
//
//					prd_f += df1 + df2 * (SUBSEARCH_RANGE - 1);
//				}
//			}
//		}
//	}
//
//	cbuf_p = cbuf;
//	j = (SEARCH_RANGE * SUBSEARCH_RANGE) >> 1;
//
//	for (i = 0; i < SEARCH_RANGE * SUBSEARCH_RANGE; i++) {
//		if (cbuf_p[i] < cbuf_p[j]) {
//			j = i;
//		}
//	}
//
//	i = (j / SUBSEARCH_RANGE) - (SEARCH_RANGE >> 1);
//	j = (j % SUBSEARCH_RANGE) - (SUBSEARCH_RANGE >> 1);
//	y = coef_p[pos1] + i;
//	x = coef_p[pos2] - i - j;
//
//	if (y < -MAX_COEF) y = -MAX_COEF;
//	else if (y > MAX_COEF) y = MAX_COEF;
//
//	if (x < -MAX_COEF) x = -MAX_COEF;
//	else if (x > MAX_COEF) x = MAX_COEF;
//
//	i = y - coef_p[pos1];
//	j = x - coef_p[pos2];
//
//	if (i != 0 || j != 0) {
//		for (y = 0; y < enc->height; y++) {
//			class_p = enc->class[y];
//
//			for (x = 0; x < enc->width; x++) {
//				if (cl == *class_p++) {
//					roff1_p = roff1_pp[y][x];
//                    roff2_p = roff2_pp[y][x];
//
//					org_p = &enc->org[1][y][x];
//					enc->prd[y][x] += org_p[roff1_p[pos1]] * i + org_p[roff2_p[pos2]] * j;
//				}
//			}
//		}
//
//		coef_p[pos1] += i;
//		coef_p[pos2] += j;
//	}
}

/*------------------------- optimize_predictor ----------------------*
 |  Function optimize_predictor
 |
 |  Purpose: Two coefficients are randomly chosen
 |			 and optimized together
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN/OUT)
 |
 |  Returns:  cost_t	--> Returns the new cost value
 *-------------------------------------------------------------------*/
cost_t optimize_predictor(ENCODER *enc) {
    int cl, k, pos1, pos2;
#ifndef RAND_MAX
#  define RAND_MAX 32767
#endif

    for (cl = 0; cl < enc->num_class; cl++) {
        for (k = 0; k < enc->full_prd_order; k++) {
            retry:
            pos1 = (int) (((double) random() * enc->full_prd_order) / (RAND_MAX + 1.0));
            pos2 = (int) (((double) random() * enc->full_prd_order) / (RAND_MAX + 1.0));

            if (pos1 == pos2) goto retry;

            optimize_coef(enc, cl, pos1, pos2);
        }
    }

    predict_region(enc, 0, 0, 0, 0, enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH]);
    return (calc_cost(enc, 0, 0, 0, 0, enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH]));
}

int putbits(FILE *fp, int n, uint x) {
    static int bitpos = 8;
    static uint bitbuf = 0;
    int bits;

    bits = n;

    if (bits <= 0) return (0);

    while (n >= bitpos) {
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

/*-------------------------- remove_emptyclass ----------------------*
 |  Function remove_emptyclass
 |
 |  Purpose: Removes the empty classes and updates
 |			 the dependent variables
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN/OUT)
 |
 |  Returns:  void
 *-------------------------------------------------------------------*/
void remove_emptyclass(ENCODER *enc) {
    int cl, i, k, v, u, t, s;

    for (cl = 0; cl < enc->num_class; cl++) {
        enc->mtfbuf[cl] = 0;
    }

    for (v = 0; v < enc->vu[HEIGHT]; v += MIN_BSIZE) {
        for (u = 0; u < enc->vu[WIDTH]; u += MIN_BSIZE) {
            for (t = 0; t < enc->ts[HEIGHT]; t += MIN_BSIZE) {
                for (s = 0; s < enc->ts[WIDTH]; s += MIN_BSIZE) {
                    cl = enc->class[v][u][t][s];
                    enc->mtfbuf[cl]++;
                }
            }
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

    for (v = 0; v < enc->vu[HEIGHT]; v++) {
        for (u = 0; u < enc->vu[WIDTH]; u++) {
            for (t = 0; t < enc->ts[HEIGHT]; t++) {
                for (s = 0; s < enc->ts[WIDTH]; s++) {
                    i = enc->class[v][u][t][s];
                    enc->class[v][u][t][s] = (char) enc->mtfbuf[i];
                }
            }
        }
    }

    for (i = cl = 0; i < enc->num_class; i++) {
        if (enc->mtfbuf[i] < 0) continue;

        if (cl != i) {
            for (k = 0; k < enc->full_prd_order; k++) {
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

/*---------------------------- write_header ---------------------------*
 |  Function write_header
 |
 |  Purpose: Writes the header of the file
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN/OUT)
 |		prd_order		--> Variable that store the intra prediction order (IN)
 |		mi_prd_order    --> Array that stores the LF prediction orders (IN)
 |		fp				--> File to write into (IN)
 |
 |  Returns:  int		--> Returns the number of bits used
 *---------------------------------------------------------------------*/
int write_header(ENCODER *enc, int prd_order, int mi_prd_order[4], FILE *fp) {
    int bits, i;

    bits =  putbits(fp, 16, MAGIC_NUMBER);
    bits += putbits(fp, 16, (uint) enc->vu[HEIGHT]);
    bits += putbits(fp, 16, (uint) enc->vu[WIDTH]);
    bits += putbits(fp, 16, (uint) enc->ts[HEIGHT]);
    bits += putbits(fp, 16, (uint) enc->ts[WIDTH]);
    bits += putbits(fp, 16, (uint) enc->maxval);
    bits += putbits(fp, 6,  (uint) enc->depth);
    bits += putbits(fp, 4,  1);	/* number of components (1 = monochrome) */
    bits += putbits(fp, 6,  (uint) enc->num_group);

    bits += putbits(fp, 8, (uint) prd_order);

    for (i = 0; i < 4; i++) {
        bits += putbits(fp, 8, (uint) mi_prd_order[i]);
    }

    bits += putbits(fp, 8,  (uint) enc->delta);
    bits += putbits(fp, 6,  (uint) enc->num_pmodel - 1);
    bits += putbits(fp, 4,  (uint) enc->coef_precision - 1);
    bits += putbits(fp, 4,  (uint) enc->pm_accuracy + 1);
    bits += putbits(fp, 1,  (uint) enc->f_huffman);
    bits += putbits(fp, 1,  (enc->quadtree_depth < 0)? 0 : 1);

    return (bits);
}

/*---------------------------- write_class ----------------------------*
 |  Function write_class
 |
 |  Purpose: Writes the number of classes for the current frame
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN/OUT)
 |		fp				--> File to write into (IN)
 |
 |  Returns:  int		--> Returns the number of bits used
 *---------------------------------------------------------------------*/
int write_class(ENCODER *enc, FILE *fp) {
    int bits;

    bits = putbits(fp, 8, (uint) enc->num_class);
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
    putbits(fp, m, (uint) v);

    return (bits + m);
}

void set_qtindex(ENCODER *enc, int *index, uint *hist, int *numidx, int tlv, int tlu, int tlt, int tls, int blksize,
                 int width, int level) {
    int i, cl, v, u, t, s, ctx, last_width;
    char ****qtmap;

    if (tlv >= enc->vu[HEIGHT] || tlu >= enc->vu[WIDTH] || tlt >= enc->ts[HEIGHT] || tls >= enc->ts[WIDTH]) return;

    if (level > 0) {
        /* context modeling for quad-tree flag */
        ctx = 0;
        qtmap = enc->qtmap[level - 1];

        v = (tlv / MIN_BSIZE) >> level;
        u = (tlu / MIN_BSIZE) >> level;
        t = (tlt / MIN_BSIZE) >> level;
        s = (tls / MIN_BSIZE) >> level;

        if (v > 0 && qtmap[v - 1][u][t][s] == 1) ctx++;

        if (u > 0 && qtmap[v][u - 1][t][s] == 1) ctx++;

        if (t > 0) {
            if (qtmap[v][u][t - 1][s] == 1) ctx++;
            if (tls + blksize < width && qtmap[v][u][t - 1][s + 1] == 1) ctx++;
        }

        if (s > 0 && qtmap[v][u][t][s - 1] == 1) ctx++;

        ctx = ((level - 1) * 4 + ctx) << 1;

        if (qtmap[v][u][t][s] == 1) {
            ctx++;
            index[(*numidx)++] = -(ctx + 1);
            enc->qtctx[ctx]++;
            blksize >>= 1;

            // v and u
            set_qtindex(enc, index, hist, numidx, tlv, tlu, tlt, tls, blksize, width, level - 1);
            set_qtindex(enc, index, hist, numidx, tlv, tlu, tlt, tls + blksize, blksize, width, level - 1);
            set_qtindex(enc, index, hist, numidx, tlv, tlu, tlt + blksize, tls, blksize, width, level - 1);

            last_width = tls + blksize * 2;
            if (last_width >= enc->ts[WIDTH]) last_width = enc->ts[WIDTH];

            set_qtindex(enc, index, hist, numidx, tlv, tlu, tlt + blksize, tls + blksize, blksize, last_width, level - 1);

            // v and u + blksize
            set_qtindex(enc, index, hist, numidx, tlv, tlu + blksize, tlt, tls, blksize, width, level - 1);
            set_qtindex(enc, index, hist, numidx, tlv, tlu + blksize, tlt, tls + blksize, blksize, width, level - 1);
            set_qtindex(enc, index, hist, numidx, tlv, tlu + blksize, tlt + blksize, tls, blksize, width, level - 1);

            last_width = tls + blksize * 2;
            if (last_width >= enc->ts[WIDTH]) last_width = enc->ts[WIDTH];

            set_qtindex(enc, index, hist, numidx, tlv, tlu + blksize, tlt + blksize, tls + blksize, blksize, last_width, level - 1);

            // v + blksize and u
            set_qtindex(enc, index, hist, numidx, tlv + blksize, tlu, tlt, tls, blksize, width, level - 1);
            set_qtindex(enc, index, hist, numidx, tlv + blksize, tlu, tlt, tls + blksize, blksize, width, level - 1);
            set_qtindex(enc, index, hist, numidx, tlv + blksize, tlu, tlt + blksize, tls, blksize, width, level - 1);

            last_width = tls + blksize * 2;
            if (last_width >= enc->ts[WIDTH]) last_width = enc->ts[WIDTH];

            set_qtindex(enc, index, hist, numidx, tlv + blksize, tlu, tlt + blksize, tls + blksize, blksize, last_width, level - 1);

            // v + blksize and u + blksize
            set_qtindex(enc, index, hist, numidx, tlv + blksize, tlu + blksize, tlt, tls, blksize, width, level - 1);
            set_qtindex(enc, index, hist, numidx, tlv + blksize, tlu + blksize, tlt, tls + blksize, blksize, width, level - 1);
            set_qtindex(enc, index, hist, numidx, tlv + blksize, tlu + blksize, tlt + blksize, tls, blksize, width, level - 1);

            last_width = tls + blksize * 2;
            if (last_width >= enc->ts[WIDTH]) last_width = enc->ts[WIDTH];

            set_qtindex(enc, index, hist, numidx, tlv + blksize, tlu + blksize, tlt + blksize, tls + blksize, blksize, last_width, level - 1);
            return;
        }
        else {
            index[(*numidx)++] = -(ctx + 1);
            enc->qtctx[ctx]++;
        }
    }

    cl = enc->class[tlv][tlu][tlt][tls];
    mtf_classlabel(enc->class, enc->mtfbuf, tlv, tlu, tlt, tls, blksize, width, enc->num_class);
    i = enc->mtfbuf[cl];
    index[(*numidx)++] = i;
    hist[i]++;
}

/*--------------------------- encode_class ----------------------------*
 |  Function encode_class
 |
 |  Purpose: Arithmetic encode of the classes to file
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN/OUT)
 |		fp				--> File to write into (IN)
 |
 |  Returns:  int		--> Returns the number of bits used
 *---------------------------------------------------------------------*/
int encode_class(FILE *fp, ENCODER *enc) {
    int i, j, k, numidx, blksize, level, v, u, t, s, ctx, bits, *index;
    uint *hist;
    cost_t cost;

#ifndef OPT_SIDEINFO
    if (fp == NULL) return(0);
#endif
    if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
        level = enc->quadtree_depth;
        blksize = MAX_BSIZE;
        numidx = 0;

        v = enc->vu[HEIGHT] / MIN_BSIZE;
        u = enc->vu[WIDTH] / MIN_BSIZE;
        t = enc->ts[HEIGHT] / MIN_BSIZE;
        s = enc->ts[WIDTH] / MIN_BSIZE;

        for (k = 0; k <= level; k++) {
            numidx += v * u * t * s;

            v = (v + 1) >> 1;
            u = (u + 1) >> 1;
            t = (t + 1) >> 1;
            s = (s + 1) >> 1;
        }

        for (k = 0; k < QUADTREE_DEPTH << 3; k++) {
            enc->qtctx[k] = 0;
        }
    }
    else {
        level = 0;
        blksize = BASE_BSIZE;
        numidx = enc->vu[HEIGHT] * enc->vu[WIDTH] * enc->ts[HEIGHT] * enc->ts[WIDTH] /
                 (BASE_BSIZE * BASE_BSIZE * BASE_BSIZE * BASE_BSIZE);
    }

    hist = (uint *) alloc_mem(enc->num_class * sizeof(uint));
    index = (int *) alloc_mem(numidx * sizeof(int));

    for (i = 0; i < enc->num_class; i++) {
        hist[i] = 0;
        enc->mtfbuf[i] = i;
    }

    numidx = 0;

    for (v = 0; v < enc->vu[HEIGHT]; v += blksize) {
        for (u = 0; u < enc->vu[WIDTH]; u += blksize) {
            for (t = 0; t < enc->ts[HEIGHT]; t += blksize) {
                for (s = 0; s < enc->ts[WIDTH]; s += blksize) {
                    set_qtindex(enc, index, hist, &numidx, v, u, t, s, blksize, enc->ts[WIDTH], level);
                }
            }
        }
    }

    bits = 0;

    if (enc->f_huffman == 1) {    /* Huffman */
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
                    bits += (int) (enc->class_cost[i]);
                }
            }
        }
        else {    /* actually encode */
            for (i = 0; i < enc->num_class; i++) {
                bits += putbits(fp, 4, (uint) (vlc->len[i] - 1));
            }

            for (j = 0; j < numidx; j++) {
                i = index[j];

                if (i < 0) {
                    ctx = -(i + 1);
                    putbits(fp, 1, (uint) (ctx & 1));
                }
                else {
                    bits += putbits(fp, vlc->len[i], vlc->code[i]);
                }
            }
        }

        free_vlc(vlc);
    }
    else {            /* Arithmetic */
        PMODEL *pm;
        double p, c;
        int qtree_code[QUADTREE_DEPTH << 2], mtf_code[MAX_CLASS];

        /* context modeling for quad-tree flag */
        if (level > 0) {
            for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
                cost = INT_MAX;

                for (i = k = 0; i < 7; i++) {
                    p = qtree_prob[i];
                    c = -log(p) * (cost_t) enc->qtctx[(ctx << 1) + 1] - log(1.0 - p) * (cost_t) enc->qtctx[ctx << 1];

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
            c += (double) hist[i];
        }

        for (i = 0; i < enc->num_class; i++) {
            p = (double) hist[i] / c;

            if (p > 0.0) {
                mtf_code[i] = (int) (-log(p) / log(2.0) * ((float) PMCLASS_LEVEL / PMCLASS_MAX));

                if (mtf_code[i] >= PMCLASS_LEVEL) {
                    mtf_code[i] = PMCLASS_LEVEL - 1;
                }
            }
            else {
                mtf_code[i] = PMCLASS_LEVEL - 1;
            }

            p = exp(-log(2.0) * ((double) mtf_code[i] + 0.5) * PMCLASS_MAX / PMCLASS_LEVEL);
            enc->class_cost[i] = -log(p) / log(2.0);
            hist[i] = (uint) (p * (1 << 10));

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

            bits = (int) cost;
        }
        else {    /* actually encode */
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
                    pm->freq[(ctx << 1) + 1] = (uint) (p * (1 << 10));
                    p = 1.0 - p;
                    pm->freq[(ctx << 1)] = (uint) (p * (1 << 10));
                }

                for (i = 0; i < pm->size; i++) {
                    pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
                }
            }

            cpm->size = enc->num_class;
            cpm->freq = (uint *) alloc_mem((cpm->size * 2 + 1) * sizeof(uint));
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
                    rc_encode(fp, enc->rc, pm->cumfreq[i] - pm->cumfreq[ctx], pm->freq[i],
                              pm->cumfreq[ctx + 2] - pm->cumfreq[ctx]);
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

/*-------------------------- encode_predictor --------------------------*
 |  Function encode_predictor
 |
 |  Purpose: Arithmetic encode of the predictors to file
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN/OUT)
 |		fp				--> File to write into (IN)
 |
 |  Returns:  int		--> Returns the number of bits used
 *---------------------------------------------------------------------*/
int encode_predictor(FILE *fp, ENCODER *enc) {
    int cl, coef, sgn, k, m, min_m, bits = 0;
    cost_t cost, min_cost, t_cost;

#ifndef OPT_SIDEINFO
    if (fp == NULL) return(0);
#endif
    t_cost = 0.0;

    for (k = 0; k < enc->full_prd_order; k++) {
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

    bits = (int) t_cost;

    if (fp != NULL) {
        bits = 0;

        if (enc->f_huffman == 1) {	/* Huffman */
            for (k = 0; k < enc->full_prd_order; k++) {
                bits += putbits(fp, 4, (uint) enc->coef_m[k]);

                for (cl = 0; cl < enc->num_class; cl++) {
                    coef = enc->predictor[cl][k];
                    sgn = (coef < 0) ? 1 : 0;

                    if (coef < 0) coef = -coef;

                    bits += encode_golomb(fp, enc->coef_m[k], coef);

                    if (coef != 0) {
                        bits += putbits(fp, 1, (uint) sgn);
                    }
                }
            }
        }
        else {			/* Arithmetic */
            PMODEL *pm;
            pm = &enc->spm;

            for (k = 0; k < enc->full_prd_order; k++) {
                set_spmodel(pm, MAX_COEF + 1, enc->coef_m[k]);
                rc_encode(fp, enc->rc, (uint) enc->coef_m[k], 1, 16);

                for (cl = 0; cl < enc->num_class; cl++) {
                    coef = enc->predictor[cl][k];
                    sgn = (coef < 0) ? 1 : 0;

                    if (coef < 0) coef = -coef;

                    rc_encode(fp, enc->rc, pm->cumfreq[coef], pm->freq[coef], pm->cumfreq[pm->size]);

                    if (coef > 0) {
                        rc_encode(fp, enc->rc, (uint) sgn, 1, 2);
                    }
                }
            }

            bits = (int) enc->rc->code;
            enc->rc->code = 0;
        }
    }

    return (bits);
}

/*-------------------------- encode_threshold --------------------------*
 |  Function encode_threshold
 |
 |  Purpose: Arithmetic encode of the thresholds to file
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN/OUT)
 |		fp				--> File to write into (IN)
 |
 |  Returns:  int		--> Returns the number of bits used
 *---------------------------------------------------------------------*/
int encode_threshold(FILE *fp, ENCODER *enc) {
    int cl, gr, i, k, m, min_m, bits = 0;
    cost_t cost, min_cost;
    PMODEL *pm;

#ifndef OPT_SIDEINFO
    if (fp == NULL) return(0);
#endif
    if (enc->f_huffman == 1) {    /* Huffman */
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
                enc->th_cost[i] = ((i - 1) >> min_m) + min_m + 1 + 1;
            }

            bits = (int) min_cost;
        }
        else {
            bits = putbits(fp, 4, (uint) min_m);

            for (cl = 0; cl < enc->num_class; cl++) {
                k = 0;

                for (gr = 1; gr < enc->num_group; gr++) {
                    i = enc->th[cl][gr - 1] - k;

                    if (i == 0) {
                        bits += putbits(fp, 1, 0);
                    }
                    else {
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
                    bits += putbits(fp, k, (uint) pm->id);
                }
            }
        }
    }
    else {            /* Arithmetic */
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
                    p = (double) pm->freq[i] / (pm->cumfreq[pm->size - k]);
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

            bits = (int) min_cost;
        }
        else {
            rc_encode(fp, enc->rc, (uint) min_m, 1, 16);

            for (cl = 0; cl < enc->num_class; cl++) {
                k = 0;

                for (gr = 1; gr < enc->num_group; gr++) {
                    i = enc->th[cl][gr - 1] - k;
                    rc_encode(fp, enc->rc, pm->cumfreq[i], pm->freq[i], pm->cumfreq[pm->size - k]);
                    k += i;

                    if (k > MAX_UPARA) break;
                }
            }

            if (enc->num_pmodel > 1) {
                for (gr = 0; gr < enc->num_group; gr++) {
                    pm = enc->pmlist[gr];
                    rc_encode(fp, enc->rc, (uint) pm->id, 1, (uint) enc->num_pmodel);
                }
            }

            bits = (int) enc->rc->code;
            enc->rc->code = 0;
        }
    }

    return (bits);
}

/*---------------------------- encode_image ---------------------------*
 |  Function encode_image
 |
 |  Purpose: Arithmetic encode of the residue to file
 |
 |  Parameters:
 |      enc		 		--> Encoder struct (IN/OUT)
 |		fp				--> File to write into (IN)
 |
 |  Returns:  int		--> Returns the number of bits used
 *---------------------------------------------------------------------*/
int encode_image(FILE *fp, ENCODER *enc) {
    int v, u, t, s, e, prd, base, bits, gr, cumbase;
    PMODEL *pm;

    bits = 0;
    if (enc->f_huffman == 1) {    /* Huffman */
        VLC *vlc;

        for (v = 0; v < enc->vu[HEIGHT]; v++) {
            for (u = 0; u < enc->vu[WIDTH]; u++) {
                for (t = 0; t < enc->ts[HEIGHT]; t++) {
                    for (s = 0; s < enc->ts[WIDTH]; s++) {
                        gr = enc->group[v][u][t][s];
                        e = enc->encval[v][u][t][s];
                        pm = enc->pmlist[gr];
                        vlc = &enc->vlcs[gr][pm->id];
                        bits += putbits(fp, vlc->len[e], vlc->code[e]);
                    }
                }
            }
        }
    }
    else {            /* Arithmetic */
        for (v = 0; v < enc->vu[HEIGHT]; v++) {
            for (u = 0; u < enc->vu[WIDTH]; u++) {
                for (t = 0; t < enc->ts[HEIGHT]; t++) {
                    for (s = 0; s < enc->ts[WIDTH]; s++) {
                        gr = enc->group[v][u][t][s];
                        prd = enc->prd[v][u][t][s];

                        if (prd < 0) prd = 0;
                        else if (prd > enc->maxprd) prd = enc->maxprd;

                        e = enc->encval[v][u][t][s];
                        base = enc->bconv[prd];
                        pm = enc->pmlist[gr] + enc->fconv[prd];
                        cumbase = pm->cumfreq[base];

                        rc_encode(fp, enc->rc, pm->cumfreq[base + e] - cumbase, pm->freq[base + e],
                                  pm->cumfreq[base + enc->maxval + 1] - cumbase);
                    }
                }
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

void print_results(FILE *res, const int *vu, const int *ts, int header, int class_info, int predictors, int thresholds,
                   int errors) {
    int rate = 0, total_bits = 0;

    printf("\n---------------------------------\n");
    printf("Header info.\t :%10d bits\n", header);
    fprintf(res, "Header info.\t :%10d bits\n", header);
    fprintf(res, "Bits\t\tBpp\n");

    // Results output
    printf("---------------------------------\n");
    printf("class info\t :%10d bits\n", class_info);
    printf("predictors\t :%10d bits\n", predictors);
    printf("thresholds\t :%10d bits\n", thresholds);
    printf("pred. errors\t :%10d bits\n", errors);

    rate = class_info + predictors + thresholds + errors;

    total_bits += rate;

    printf("---------------------------------\n");
    printf("total LF:%10d bits\n", rate);
    printf("LF coding rate:%10.5f bpp\n", (double) rate / (vu[HEIGHT] * vu[WIDTH] * ts[HEIGHT] * ts[WIDTH]));
    fprintf(res, "%10d\t%10.5f\n", rate, (double) rate / (vu[HEIGHT] * vu[WIDTH] * ts[HEIGHT] * ts[WIDTH]));

    total_bits += header;

    printf("---------------------------------\n");
    printf("total\t\t :%10d bits\n", total_bits);
    printf("total coding rate:%10.5f bpp\n", (double) total_bits / (vu[HEIGHT] * vu[WIDTH] * ts[HEIGHT] * ts[WIDTH]));
}

/*--------------------------- total_variation ---------------------------*
 |  Function total_variation
 |
 |  Purpose: Calculates the total variation of an image
 |
 |  Parameters:
 |		img				--> Image to analyse (IN)
 |
 |  Returns:  double	--> Returns the total variation
 *----------------------------------------------------------------------*/
double total_variation(LF4D *lf) {
    int v, u, t, s;

    double tv = 0.0, dt, ds;

    for (v = 0; v < lf->v; v++) {
        for (u = 0; u < lf->u; u++) {
            for (t = 0; t < lf->t; t++) {
                for (s = 0; s < lf->s; s++) {
                    dt = (t == 0) ? 0 : (lf->val[v][u][t][s] - lf->val[v][u][t - 1][s]) * (lf->val[v][u][t][s] - lf->val[v][u][t - 1][s]);
                    ds = (s == 0) ? 0 : (lf->val[v][u][t][s] - lf->val[v][u][t][s - 1]) * (lf->val[v][u][t][s] - lf->val[v][u][t][s - 1]);

                    tv = tv + sqrt(dt + ds);
                }
            }
        }
    }

    return tv;
}

/*--------------------------- sparseness_index ----------------------------*
 |  Function sparseness_index
 |
 |  Purpose: Calculates the percentage of used values in the dynamic range
 |
 |  Parameters:
 |		infile		--> Image/Sequence to analyse (IN)
 |		height		--> Height of the video (IN)
 |		width		--> Width of the video (IN)
 |		frame		--> Frame of the video to copy (IN)
 |		depth		--> Image dynamic range (bpp) (IN)
 |		endianness	--> YUV endianness, for depth > 8 bpp (IN)
 |		chroma  	--> YUV chroma type (IN)
 |
 |  Returns: double	--> Returns the sparseness index
 *-------------------------------------------------------------------------*/
double sparseness_index(char *infile, int *mi_size, int *vu, int depth, int endianness, int chroma, int format) {
    int k, v, u, t, s, L = 0;
    unsigned long int *hist = (unsigned long int *) alloc_mem(sizeof(unsigned long int) * ((unsigned long int) pow(2, depth)));
    int max = 0, min = ((int) pow(2, depth) - 1);
    char *extra_info = NULL;
    LF4D *lf = NULL;

    for (k = 0; k < (int) pow(2, depth); k++) {
        hist[k] = 0;
    }

    // Operations to obtain the lookup table
    lf = read_yuv(infile, vu[HEIGHT], vu[WIDTH], mi_size[HEIGHT], mi_size[WIDTH], depth, endianness, chroma, format);

    for (v = 0; v < lf->v; v++) {
        for (u = 0; u < lf->u; u++) {
            for (t = 0; t < lf->t; t++) {
                for (s = 0; s < lf->s; s++) {
                    hist[lf->val[v][u][t][s]]++;
                }
            }
        }
    }

    safefree_lf4d(&lf);

    for (k = 0; k < (int) pow(2, depth); k++) {
        if (hist[k] != 0) {
            L++;

            if (k < min) min = k;
            if (k > max) max = k;
        }
    }

    free(hist);
    safefree((void **) &extra_info);

    return (1 - (L / (double) (1 + max - min))) * 100.0;
}

/*--------------------------- histogram_check ---------------------------*
 |  Function histogram_check
 |
 |  Purpose: Checks the image histogram and decide if histogram packing is used
 |
 |  Parameters:
 |		infile		--> Image/Sequence to analyse (IN)
 |		height		--> Height of the video (IN)
 |		width		--> Width of the video (IN)
 |		frame		--> Frame of the video to copy (IN)
 |		depth		--> Image dynamic range (bpp) (IN)
 |		endianness	--> YUV endianness, for depth > 8 bpp (IN)
 |		chroma  	--> YUV chroma type (IN)
 |
 |  Returns:  int*	--> Returns the lookup table
 *----------------------------------------------------------------------*/
int* histogram_check(char *infile, int *mi_size, int *vu, int depth, int endianness, int chroma, int format) {
    int k, t, s, u, v;
    int used_values = 0;
    unsigned long int *hist = (unsigned long int *) alloc_mem(sizeof(unsigned long int) * ((unsigned long int) pow(2, depth)));
    int *forward_table = (int *) alloc_mem(sizeof(int) * (int) (pow(2, depth)));
    int *table = (int *) alloc_mem(sizeof(int) * (int) pow(2, depth));
    double tv_original = 0.0, tv_packed = 0.0;
    char *extra_info = NULL;
    LF4D *aux_lf, *lf = NULL;

    for (k = 0; k < (int) pow(2, depth); k++) {
        hist[k] = 0;
        forward_table[k] = 0;
    }

    // Operations to obtain the lookup table
    lf = read_yuv(infile, vu[HEIGHT], vu[WIDTH], mi_size[HEIGHT], mi_size[WIDTH], depth, endianness, chroma, format);

    for (v = 0; v < lf->v; v++) {
        for (u = 0; u < lf->u; u++) {
            for (t = 0; t < lf->t; t++) {
                for (s = 0; s < lf->s; s++) {
                    hist[lf->val[v][u][t][s]]++;
                }
            }
        }
    }

    safefree_lf4d(&lf);

    for (k = 0; k < (int) pow(2, depth); k++) {
        if (hist[k] != 0) {
            forward_table[k] = 1;
        }
    }

    free(hist);

    aux_lf = alloc_lf4d(vu[HEIGHT], vu[WIDTH], mi_size[HEIGHT], mi_size[WIDTH], (int) (pow(2, depth) - 1));

    // Produces the actual forward table
    for (k = 0; k < (int) pow(2, depth); k++) {
        if (forward_table[k] != 0) {
            table[k] = used_values;
            used_values++;
        }
        else {
            table[k] = -1;
        }
    }

    // Perform the histogram packing of the image
    lf = read_yuv(infile, vu[HEIGHT], vu[WIDTH], mi_size[HEIGHT], mi_size[WIDTH], depth, endianness, chroma, format);

    for (v = 0; v < lf->v; v++) {
        for (u = 0; u < lf->u; u++) {
            for (t = 0; t < lf->t; t++) {
                for (s = 0; s < lf->s; s++) {
                    aux_lf->val[v][u][t][s] = (img_t) table[lf->val[v][u][t][s]];
                }
            }
        }
    }

    tv_original += total_variation(lf);
    tv_packed += total_variation(aux_lf);

    safefree_lf4d(&lf);
    safefree_lf4d(&aux_lf);

    // Check if the histogram packing lowers the total variation
    if (tv_original <= tv_packed) {
        safefree((void **) &forward_table);
    }

    free(table);
    safefree((void **) &extra_info);

    return forward_table;
}

/*------------------------- histogram_packing -------------------------*
 |  Function histogram_packing
 |
 |  Purpose: Peforms the histogram packing and checks if the
 |			 total variation is lower
 |
 |  Parameters:
 |		img				--> Image to pack (IN/OUT)
 |		forward_table	--> Lookup table to use (IN)
 |
 |  Returns:  void
 *----------------------------------------------------------------------*/
void histogram_packing(LF4D *lf, const int *forward_table, int depth) {
    int k, v, u, t, s;
    int used_values = 0;
    int *table = (int *) alloc_mem(sizeof(int) * (int) pow(2, depth));

    // Produces the actual forward table
    for (k = 0; k < (int) pow(2, depth); k++) {
        if (forward_table[k] != 0) {
            table[k] = used_values;
            used_values++;
        }
        else {
            table[k] = -1;
        }
    }

    // Perform the histogram packing of the image
    for (v = 0; v < lf->v; v++) {
        for (u = 0; u < lf->u; u++) {
            for (t = 0; t < lf->t; t++) {
                for (s = 0; s < lf->s; s++) {
                    lf->val[v][u][t][s] = (img_t) table[lf->val[v][u][t][s]];
                }
            }
        }
    }

    lf->maxval = used_values - 1;

    free(table);
}

// inserts into subject[] at position pos
char *append(char *subject, char *insert, int pos) {
    uint i = 0;

    char *buf = (char *) alloc_mem(sizeof(char) * (strlen(subject) + strlen(insert) + 1));
    buf[0] = '\0';

    for (i = 0; i < ((uint) pos < strlen(subject) ? (uint) pos : strlen(subject)); i++) {
        buf[i] = subject[i];
    }
    buf[i] = '\0';

    int len = (int) strlen(buf);
    strcpy(buf + len, insert); // copy all of insert[] at the end

    len += strlen(insert);  // increase the length by length of insert[]
    strcpy(buf + len, subject + pos); // copy the rest

    free(subject);

    return(buf);
}

/*------------------------- encode_lookuptable -------------------------*
 |  Function encode_lookuptable
 |
 |  Purpose: Encodes the lookup table using RLE
 |
 |  Parameters:
 |		fp				--> File to write into (IN/OUT)
 |		forward_table	--> Lookup table to encode (IN)
 |		size			--> Length of the lookup table (IN)
 |
 |  Returns:  int		--> Number of written bits
 *----------------------------------------------------------------------*/
int encode_lookuptable(FILE *fp, const int *forward_table, int size) {
    // Statistics
    int y = 0;
    int count = 0;
    int bit = 0;
    int bits;
    char *new_str = (char *) alloc_mem(sizeof(char));
    new_str[0] = '\0';

    for (y = 0; y < size; y++) {
        if (count == 0) {
            bit = forward_table[y];
            count++;
        }
        else {
            if (bit == forward_table[y]) {
                count++;
            }
            else {
                if (bit == 0) {
                    if (count <= 126) {
                        new_str = cat_str(new_str, cat_str("0", int2bin(count - 1, 7), 0), 1);
                    }
                    else if (count <= 382) {
                        new_str = cat_str(new_str, cat_str("01111110", int2bin(count - 127, 8), 0), 1);
                    }
                    else {
                        new_str = cat_str(new_str, cat_str("01111111", int2bin(count - 383, 16), 0), 1);
                    }
                }
                else {
                    if (count <= 126) {
                        new_str = cat_str(new_str, cat_str("1", int2bin(count - 1, 7), 0), 1);
                    }
                    else if (count <= 383) {
                        new_str = cat_str(new_str, cat_str("11111110", int2bin(count - 127, 8), 0), 1);
                    }
                    else {
                        new_str = cat_str(new_str, cat_str("11111111", int2bin(count - 383, 16), 0), 1);
                    }
                }
                count = 0;
            }
        }
    }

    if (count != 0) {
        if (bit == 0) {
            if (count <= 126) {
                new_str = cat_str(new_str, cat_str("0", int2bin(count - 1, 7), 0), 1);
            }
            else if (count <= 383) {
                new_str = cat_str(new_str, cat_str("01111110", int2bin(count - 127, 8), 0), 1);
            }
            else {
                new_str = cat_str(new_str, cat_str("01111111", int2bin(count - 383, 16), 0), 1);
            }
        }
        else {
            if (count <= 126) {
                new_str = cat_str(new_str, cat_str("1", int2bin(count - 1, 7), 0), 1);
            }
            else if (count <= 383) {
                new_str = cat_str(new_str, cat_str("11111110", int2bin(count - 127, 8), 0), 1);
            }
            else {
                new_str = cat_str(new_str, cat_str("11111111", int2bin(count - 383, 16), 0), 1);
            }
        }
    }

    // Number of bytes to write in the file
    int bytes = (int) (strlen(new_str) / 8.0);

    // Separates the bytes in the string
    for (y = 1; y < bytes; y++) {
        new_str = append(new_str, " ", 8 * y + y - 1);
    }

    char* sEnd = new_str;

    bits = putbits(fp, 16, (uint) bytes);

    // Write to file
    for (y = 0; y < bytes; y++) {
        bits += putbits(fp, 8, (uint) strtol(sEnd, &sEnd, 2));
    }

    safefree((void **) &new_str);

    return bits;
}

// Write classes to file
void debug_predictors(ENCODER *enc) {
    if (enc->debug_path != NULL) {
        char *predictor_file = (char *) alloc_mem(
                (strlen(enc->debug_path) + strlen("/predictors.txt\0") + 1) * sizeof(char));

        sprintf(predictor_file, "%s/predictors.txt", enc->debug_path);

        FILE *debug_fp = fopen(predictor_file, "w");

        fprintf(debug_fp, "#CLASS\t\tCOEFFICIENTS");

        for (int cl = 0; cl < enc->num_class; cl++) {
            fprintf(debug_fp, "\n%2d\t\t\t", cl);
            for (int p = 0; p < enc->full_prd_order; p++) {
                fprintf(debug_fp, " %3d", enc->predictor[cl][p]);
            }
        }

        fclose(debug_fp);
        free(predictor_file);
    }
}

// Write partition to file
void debug_partition(ENCODER *enc, int endianness) {
    int d, k, i, v, u, t, s, a, b;
    unsigned short byte;

    img_t *qt_lf_ptr[3];

    if (enc->debug_path != NULL) {
        char *partition_img = (char *) alloc_mem(
                (strlen(enc->debug_path) + strlen("/partition_0000x0000_10bpp_LE_YUV444p.yuv\0")) * sizeof(char));

        LF4D *qt_lf[3];
        for (k = 0; k < 3; k++) {
            qt_lf[k] = alloc_lf4d(enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH], enc->maxval);
            qt_lf_ptr[k] = &qt_lf[k]->val[0][0][0][0];
        }

        uint scale = (uint) floor(pow(2, enc->depth) / enc->num_class);

        for (v = 0; v < enc->vu[HEIGHT]; v++) {
            for (u = 0; u < enc->vu[WIDTH]; u++) {
                for (t = 0; t < enc->ts[HEIGHT]; t++) {
                    for (s = 0; s < enc->ts[WIDTH]; s++) {
                        *qt_lf_ptr[0]++ = enc->org[v][u][t][s];

                        *qt_lf_ptr[1]++ = enc->class[v][u][t][s] * scale;
                        *qt_lf_ptr[2]++ = enc->class[v][u][t][s] * scale;
                    }
                }
            }
        }

        if (enc->quadtree_depth > 0) {
            uint blksize = MAX_BSIZE;

            // TODO: Rever e fazer bonito
//            for (d = enc->quadtree_depth - 1; d >= 0; d--) {
//                for (v = 0; v < enc->vu[HEIGHT]; v += blksize) {
//                    for (u = 0; u < enc->vu[WIDTH]; u += blksize) {
//                        for (t = 0; t < enc->ts[HEIGHT]; t += blksize) {
//                            for (s = 0; s < enc->ts[WIDTH]; s += blksize) {
//                                for (a = v; a < (v + blksize < enc->vu[HEIGHT] ? v + blksize : enc->vu[HEIGHT]); a++) {
//                                    for (b = u; b < (u + blksize < enc->vu[WIDTH] ? u + blksize : enc->vu[WIDTH]); b++) {
//                                        if ((enc->qtmap[d][v / blksize][u / blksize][t / blksize][s / blksize] == 0 &&
//                                             d == enc->quadtree_depth - 1) ||
//                                            (enc->qtmap[d][v / blksize][u / blksize][t / blksize][s / blksize] == 0 &&
//                                             d < enc->quadtree_depth - 1 &&
//                                             enc->qtmap[d + 1][v / (blksize * 2)][u / (blksize * 2)][t / (blksize * 2)][
//                                                     s / (blksize * 2)] == 1)) {
//                                            if (t + blksize - 1 < enc->ts[HEIGHT]) {
//                                                for (i = s; i < (s + blksize < enc->ts[WIDTH] ? s + blksize
//                                                                                              : enc->ts[WIDTH]); i++) {
//                                                    qt_lf[0]->val[a][b][t + blksize - 1][i] = (img_t) enc->maxval;
//                                                    qt_lf[1]->val[a][b][t + blksize - 1][i] = (img_t) enc->maxval;
//                                                    qt_lf[2]->val[a][b][t + blksize - 1][i] = (img_t) enc->maxval;
//                                                }
//                                            }
//                                            if (s + blksize - 1 < enc->ts[WIDTH]) {
//                                                for (i = t; i < (t + blksize < enc->ts[HEIGHT] ? t + blksize
//                                                                                               : enc->ts[HEIGHT]); i++) {
//                                                    qt_lf[0]->val[a][b][t][i + blksize - 1] = (img_t) enc->maxval;
//                                                    qt_lf[1]->val[a][b][t][i + blksize - 1] = (img_t) enc->maxval;
//                                                    qt_lf[2]->val[a][b][t][i + blksize - 1] = (img_t) enc->maxval;
//                                                }
//                                            }
//                                        }
//                                        if (d == 0 &&
//                                            enc->qtmap[d][v / blksize][u / blksize][t / blksize][s / blksize] == 1) {
//                                            if (t + blksize / 2 - 1 < enc->ts[HEIGHT]) {
//                                                for (i = s; i < (s + blksize < enc->ts[WIDTH] ? s + blksize
//                                                                                              : enc->ts[WIDTH]); i++) {
//                                                    qt_lf[0]->val[a][b][t + blksize / 2 - 1][i] = (img_t) enc->maxval;
//                                                    qt_lf[1]->val[a][b][t + blksize / 2 - 1][i] = (img_t) enc->maxval;
//                                                    qt_lf[2]->val[a][b][t + blksize / 2 - 1][i] = (img_t) enc->maxval;
//                                                }
//                                            }
//                                            if (t + blksize - 1 < enc->ts[HEIGHT]) {
//                                                for (i = s; i < (s + blksize < enc->ts[WIDTH] ? s + blksize
//                                                                                              : enc->ts[WIDTH]); i++) {
//                                                    qt_lf[0]->val[a][b][t + blksize - 1][i] = (img_t) enc->maxval;
//                                                    qt_lf[1]->val[a][b][t + blksize - 1][i] = (img_t) enc->maxval;
//                                                    qt_lf[2]->val[a][b][t + blksize - 1][i] = (img_t) enc->maxval;
//                                                }
//                                            }
//                                            if (s + blksize / 2 - 1 < enc->ts[WIDTH]) {
//                                                for (i = t; i < (t + blksize < enc->ts[HEIGHT] ? t + blksize
//                                                                                               : enc->ts[HEIGHT]); i++) {
//                                                    qt_lf[0]->val[a][b][i][s + blksize / 2 - 1] = (img_t) enc->maxval;
//                                                    qt_lf[1]->val[a][b][i][s + blksize / 2 - 1] = (img_t) enc->maxval;
//                                                    qt_lf[2]->val[a][b][i][s + blksize / 2 - 1] = (img_t) enc->maxval;
//                                                }
//                                            }
//                                            if (s + blksize - 1 < enc->ts[WIDTH]) {
//                                                for (i = t; i < (t + blksize < enc->ts[HEIGHT] ? t + blksize
//                                                                                               : enc->ts[HEIGHT]); i++) {
//                                                    qt_lf[0]->val[a][b][i][t + blksize - 1] = (img_t) enc->maxval;
//                                                    qt_lf[1]->val[a][b][i][t + blksize - 1] = (img_t) enc->maxval;
//                                                    qt_lf[2]->val[a][b][i][t + blksize - 1] = (img_t) enc->maxval;
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
        }

        // Quadtree partition image name
        sprintf(partition_img, "%s/partition_%dx%d_%dbpp_LE_YUV444p.yuv", enc->debug_path,
                enc->vu[WIDTH] * enc->ts[WIDTH], enc->vu[HEIGHT] * enc->ts[HEIGHT], enc->depth);

        // Write image to file
        for (k = 0; k < 3; k++) {
            write_yuv(qt_lf[k], partition_img, enc->depth, endianness, SAI);
        }

        free(partition_img);

        for (k = 0; k < 3; k++) {
            safefree_lf4d(&qt_lf[k]);
        }

        // Quadtree map
//        if (enc->quadtree_depth > 0) {
//            char qtmap[1000];
//            sprintf(qtmap, "%s/qtmap.txt", enc->debug_path);
//
//            debug_fp = fopen(qtmap, "w");
//
//            uint x, y, xx, yy;
//
//            yy = (enc->height + MAX_BSIZE - 1) / MAX_BSIZE;
//            xx = (enc->width + MAX_BSIZE - 1) / MAX_BSIZE;
//
//            fprintf(debug_fp, "#LEVEL\t\tPARTITION\n");
//
//            for (i = enc->quadtree_depth - 1; i >= 0; i--) {
//                fprintf(debug_fp, "%d\t\t\t", i);
//                for (y = 0; y < yy; y++) {
//                    for (x = 0; x < xx; x++) {
//                        fprintf(debug_fp, " %d", enc->qtmap[i][y][x]);
//                    }
//                }
//
//                fprintf(debug_fp, "\n");
//
//                yy <<= 1u;
//                xx <<= 1u;
//            }
//
//            fclose(debug_fp);
//        }
    }
}

int main(int argc, char **argv) {
    // Variable declaration
    cost_t cost, min_cost, side_cost = 0;
    int i, j, k, v, u, t, s, cl, **prd_save, **th_save, **error = NULL;
    char ****class_save;
    LF4D *lf = NULL;
    ENCODER *enc = NULL;
    double elapse = 0.0;
    int f_mmse = 0;
    int f_optpred = 0;
    int f_huffman = 0;
    int do_histogram_packing = 1;
    int *forward_table = NULL;
    int quadtree_depth = QUADTREE_DEPTH;
    int num_class = NUM_CLASS;
    int num_group = NUM_GROUP;
    int mi_size[2] = {0, 0};
    int vu[2] = {0, 0};
    int prd_order = 0, mi_prd_order[4] = {0, 0, 0, 0};
    int coef_precision = COEF_PRECISION;
    int num_pmodel = NUM_PMODEL;
    int pm_accuracy = PM_ACCURACY;
    int max_iteration = MAX_ITERATION;
    int depth = DEPTH;
    int endianness = LITTLE_ENDIANNESS;
    int chroma = GRAY;
    char *chroma_name = NULL;
    char *infile, *outfile;
    char resfile[1000];
    int format = SAI;
    char *format_name = NULL;
    int debug = 0;
    char *debug_path = NULL;
    FILE *fp, *res;

    //Print results variables
    int header = 0, class_info, predictors, thresholds, errors;
    int delta = 1;

    cpu_time();
    setbuf(stdout, 0);
    infile = outfile = NULL;

    // Read input parameters
    for (i = 1; i < argc; i++) {
        if (argv[i][0] == '-') {
            switch(argv[i][1]) {
                case 'J':
                    mi_size[HEIGHT] = (int) strtol(argv[++i], NULL, 10);
                    mi_size[WIDTH] = (int) strtol(argv[++i], NULL, 10);

                    if (mi_size[HEIGHT] < 0 || mi_size[WIDTH] < 0) {
                        fprintf(stderr, "The dimensions of the micro images cannot be less than 1!\n");
                        exit(-3);
                    }

                    break;

                case 'K':
                    for (j = 0; j < 2; j++) {
                        vu[j] = (int) strtol(argv[++i], NULL, 10);

                        if (vu[j] <= 0) {
                            fprintf(stderr, "The views dimensions cannot be less than 0!\n");
                            exit(-1);
                        }
                    }

                    break;

                case 'M':
                    num_class = (int) strtol(argv[++i], NULL, 10);

                    if (num_class <= 0 || num_class > 63) {
                        num_class = NUM_CLASS;
                    }

                    break;

                case 'L':
                    // TODO: prd_order[0] and mi_prd_order can't be 1 for some reason...
                    prd_order = (int) strtol(argv[++i], NULL, 10);

                    if (prd_order <= -2 || prd_order == 0 || prd_order > 72) {
                        prd_order = INTRA_PRD_ORDER;
                    }

                    for (j = 0; j < 4; j++) {
                        mi_prd_order[j] = (int) strtol(argv[++i], NULL, 10);

                        if ((mi_prd_order[j] > 0 && mi_prd_order[j] < 2) || mi_prd_order[j] > 72) {
                            mi_prd_order[j] = MI_PRD_ORDER;
                        }
                    }

                    break;

                case 'P':
                    coef_precision = (int) strtol(argv[++i], NULL, 10);

                    if (coef_precision <= 0 || coef_precision > 16) {
                        coef_precision = COEF_PRECISION;
                    }

                    break;

                case 'V':
                    num_pmodel = (int) strtol(argv[++i], NULL, 10);

                    if (num_pmodel <= 0 || num_pmodel > 64) {
                        num_pmodel = NUM_PMODEL;
                    }

                    break;

                case 'A':
                    pm_accuracy = (int) strtol(argv[++i], NULL, 10);

                    if (pm_accuracy < -1 || pm_accuracy > 6) {
                        pm_accuracy = PM_ACCURACY;
                    }

                    break;

                case 'I':
                    max_iteration = (int) strtol(argv[++i], NULL, 10);

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

                case 'u':
                    do_histogram_packing = 0;

                    break;

                case 'D':
                    delta = (int) strtol(argv[++i], NULL, 10);

                    if (delta <= 0) {
                        delta = 1;
                    }

                    break;

                case 'b':
                    depth = (int) strtol(argv[++i], NULL, 10);

                    if (depth < 8) {
                        depth = DEPTH;
                    }
                    else if (depth > 16) {
                        depth = 16;
                    }

                    break;

                case 'E':
                    endianness = (int) strtol(argv[++i], NULL, 10);

                    if (endianness == 0) {
                        endianness = LITTLE_ENDIANNESS;
                    }
                    else if (endianness == 1) {
                        endianness = BIG_ENDIANNESS;
                    }
                    else if (endianness != LITTLE_ENDIANNESS && endianness != BIG_ENDIANNESS) {
                        endianness = LITTLE_ENDIANNESS;
                    }

                    break;

                case 'C':
                    chroma_name = argv[++i];

                    if (strcmp(chroma_name, "GRAY") == 0) {
                        chroma = GRAY;
                    }
                    else if (strcmp(chroma_name, "444") == 0) {
                        chroma = S444;
                    }
                    else if (strcmp(chroma_name, "422") == 0) {
                        chroma = S422;
                    }
                    else if (strcmp(chroma_name, "411") == 0) {
                        chroma = S411;
                    }
                    else if (strcmp(chroma_name, "420") == 0) {
                        chroma = S420;
                    }
                    else {
                        printf("Chroma format not recognised: %s.\n", chroma_name);
                        printf("Supported formats:\n");
                        printf("\tGRAY;\n"
                               "\t444;\n"
                               "\t422; --> Not yet implemented\n"
                               "\t411; --> Not yet implemented\n"
                               "\t420.\n");
                        exit(-2);
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
        printf(BANNER"\n", MRP_VERSION, MRP_VERSION_DATE);
        printf("usage: encmrp [options] infile outfile\n");
        printf("options:\n");
        printf("    -J 2 * num  Views dimensions (in pixels) [%c %c]\n", 'H', 'W');
        printf("    -K 2 * num  Dimensions of the array of views [%c %c]\n", 'H', 'W');
        printf("    -L 5 * num  Prediction order (1 * Intra, 4 * Inter) [%d %d %d %d %d]\n", INTRA_PRD_ORDER, mi_prd_order[UP], mi_prd_order[LEFT], mi_prd_order[LDIAG], mi_prd_order[RDIAG]);
        printf("    -b num  	Bit depth [%d]\n", depth);
        printf("    -E num  	Endianness: little-endian = 0, big-endian = 1. Default: %s\n", "little-endian");
        printf("    -C str  	Chroma format [%s]. Supported formats:\n", "GRAY");
        printf("\t\tGRAY;\n"
               "\t\t444;\n"
               "\t\t422; --> Not yet implemented\n"
               "\t\t411; --> Not yet implemented\n"
               "\t\t420.\n"
               "\t\t(Notice: Currently MRP is a Luma only encoder. Thus this step is used only to skip the Chromas.)\n");
        printf("    -D num  	Distance between views [%d]*\n", delta);
        printf("    -M num  	Number of predictors [%d]\n", num_class);
        printf("    -P num  	Precision of prediction coefficients (fractional bits) [%d]\n", coef_precision);
        printf("    -V num  	Number of probability models [%d]\n", num_pmodel);
        printf("    -A num  	Accuracy of probability models [%d]\n", pm_accuracy);
        printf("    -I num  	Maximum number of iterations [%d]\n", max_iteration);
        printf("    -m      	Use MMSE predictors\n");
        printf("    -h      	Use Huffman coding\n");
        printf("    -f      	Fixed block-size for adaptive prediction\n");
        printf("    -u      	Deactivate the histogram packing procedures\n");
        printf("    -o      	Further optimization of predictors (experimental)\n");
        printf("    -d  		Create extra debug output (coefficients, partitions, etc.)");
        printf("    -r str 		Light field file format [%s]. Supported formats:\n", "SAI");
        printf("\t\tMIA; --> Not yet implemented\n"
               "\t\tPVS;\n"
               "\t\tSAI.\n");
        printf("infile:     	Input file (must be in a raw YUV format)\n");
        printf("outfile:    	Output file\n");
        printf("\nNote: * stands for a mandatory option.\n");
        exit(0);
    }

    // Open output file
    fp = fileopen(outfile, "wb");

    // If the number of classes was not defined it is now
    k = vu[HEIGHT] * vu[WIDTH] * mi_size[HEIGHT] * mi_size[WIDTH];
    if (num_class < 0) {
        num_class = (int) (10.4E-5 * k + 13.8);

        if (num_class > MAX_CLASS) {
            num_class = MAX_CLASS;
        }
    }

    // If the prediction order was not defined it is now
    if (prd_order < 0) {
        prd_order = (int) (12.0E-5 * k + 17.2);

        for (i = 1; i < 8; i++) {
            if (prd_order < (i + 1) * (i + 1)) {
                prd_order = i * (i + 1);
                break;
            }
        }

        if (i >= 8) prd_order = 72;
    }

    // Create directory to store debug information
    if (debug == 1) {
        // Get debug path name
        int len = strlen(outfile);
        debug_path = (char *) alloc_mem((len + 1) * sizeof(char));
        strcpy(debug_path, outfile);
        strcpy(&debug_path[len - 4], "_enc");

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

    if (do_histogram_packing == 1) forward_table = histogram_check(infile, mi_size, vu, depth, endianness, chroma, format);

    printf("\nMRP-Video Encoder\n\n");
    // Print file characteristics to screen
    printf("%s (%d x %d x %d x %d x %d) -> %s\n", infile, vu[HEIGHT], vu[WIDTH], mi_size[HEIGHT], mi_size[WIDTH],depth, outfile);
    // Print coding parameters to screen
    printf("M = %d, P = %d, V = %d, A = %d, D = %d\n\n", num_class, coef_precision, num_pmodel, pm_accuracy, delta);
    // Print prediction parameters to screen
    if (forward_table != NULL) {
        printf("Sparseness Index, S = %3.1f%%\n\n", sparseness_index(infile, mi_size, vu, depth, endianness, chroma, format));
    }

    printf("Prediction order:\n\tFrame I: %d\n\n", prd_order);

    printf("LF Prediction order: %d %d %d %d ", mi_prd_order[UP], mi_prd_order[LEFT], mi_prd_order[LDIAG],
           mi_prd_order[RDIAG]);
    printf("(4DLF dimensions: %d x %d x %d x %d)\n\n", vu[HEIGHT], vu[WIDTH], mi_size[HEIGHT], mi_size[WIDTH]);

    char *aux = remove_ext(outfile, '.', '/');
    sprintf(resfile, "%s/res_%s.txt", dirname(aux), basename(aux));
    free(aux);

    res = fileopen(resfile, "w");
    fprintf(res, "MRP-Lenslet version %s (%s) encoding results\n", MRP_VERSION, MRP_VERSION_DATE);
    fprintf(res, "\tEncoded file: %s\n", infile);
    fprintf(res, "\tDimensions: %d x %d x %d x %d\n", vu[HEIGHT], vu[WIDTH], mi_size[HEIGHT], mi_size[WIDTH]);
    fprintf(res, "---------------------------------------------\n");

    //// Start of the encoding process
    // Read input file
    lf = read_yuv(infile, vu[HEIGHT], vu[WIDTH], mi_size[HEIGHT], mi_size[WIDTH], depth, endianness, chroma, format);

    // Perform histogram packing
    if (forward_table != NULL) histogram_packing(lf, forward_table, depth);

    // Creates new ENCODER structure
    enc = init_encoder(lf, num_class, num_group, prd_order, mi_prd_order, coef_precision, f_huffman, quadtree_depth,
                       num_pmodel, pm_accuracy, delta, depth, debug_path);

    // Set cost model
    set_cost_model(enc, f_mmse);

    // First block classification
    init_class(enc);

    // Auxiliary variables
    prd_save = (int **) alloc_2d_array(enc->num_class, enc->full_prd_order, sizeof(int));
    th_save = (int **) alloc_2d_array(enc->num_class, enc->num_group, sizeof(int));
    class_save = (char ****) alloc_4d_array(enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH], sizeof(char));

    /* 1st loop */
    // Loop type
    enc->optimize_loop = 1;
    min_cost = INT_MAX;

    for (i = j = 0; i < max_iteration; i++) {
        design_predictor(enc, f_mmse);
        optimize_group(enc);
        cost = optimize_class(enc);

        if (cost < min_cost) {
            min_cost = cost;
            j = i;

            for (v = 0; v < enc->vu[HEIGHT]; v++) {
                for (u = 0; u < enc->vu[WIDTH]; u++) {
                    for (t = 0; t < enc->ts[HEIGHT]; t++) {
                        for (s = 0; s < enc->ts[WIDTH]; s++) {
                            class_save[v][u][t][s] = enc->class[v][u][t][s];
                        }
                    }
                }
            }

            for (cl = 0; cl < enc->num_class; cl++) {
                for (k = 0; k < enc->full_prd_order; k++) {
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
    for (v = 0; v < enc->vu[HEIGHT]; v++) {
        for (u = 0; u < enc->vu[WIDTH]; u++) {
            for (t = 0; t < enc->ts[HEIGHT]; t++) {
                for (s = 0; s < enc->ts[WIDTH]; s++) {
                    enc->class[v][u][t][s] = class_save[v][u][t][s];
                }
            }
        }
    }

    for (cl = 0; cl < enc->num_class; cl++) {
        for (k = 0; k < enc->full_prd_order; k++) {
            enc->predictor[cl][k] = prd_save[cl][k];
        }

        for (k = 0; k < enc->num_group; k++) {
            enc->th[cl][k] = th_save[cl][k];
        }
    }

    set_cost_rate(enc);
    predict_region(enc, 0, 0, 0, 0, enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH]);
    cost = calc_cost(enc, 0, 0, 0, 0, enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH]);

    printf("1st optimization --> Cost: %d\n", (int) cost);

    /* 2nd loop */
    //Loop type
    enc->optimize_loop = 2;
    min_cost = INT_MAX;

    for (i = j = 0; i < max_iteration; i++) {
        if (f_optpred) {
            optimize_predictor(enc);
        }

        side_cost = encode_predictor(NULL, enc);
        optimize_group(enc);
        side_cost += encode_threshold(NULL, enc);
        cost = optimize_class(enc);
        side_cost += encode_class(NULL, enc);
        cost += side_cost;

        if (cost < min_cost) {
            min_cost = cost;
            j = i;

            if (f_optpred) {
                for (v = 0; v < enc->vu[HEIGHT]; v++) {
                    for (u = 0; u < enc->vu[WIDTH]; u++) {
                        for (t = 0; t < enc->ts[HEIGHT]; t++) {
                            for (s = 0; s < enc->ts[WIDTH]; s++) {
                                class_save[v][u][t][s] = enc->class[v][u][t][s];
                            }
                        }
                    }
                }

                for (cl = 0; cl < enc->num_class; cl++) {
                    for (k = 0; k < enc->full_prd_order; k++) {
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
        for (v = 0; v < enc->vu[HEIGHT]; v++) {
            for (u = 0; u < enc->vu[WIDTH]; u++) {
                for (t = 0; t < enc->ts[HEIGHT]; t++) {
                    for (s = 0; s < enc->ts[WIDTH]; s++) {
                        enc->class[v][u][t][s] = class_save[v][u][t][s];
                    }
                }
            }
        }

        for (cl = 0; cl < enc->num_class; cl++) {
            for (k = 0; k < enc->full_prd_order; k++) {
                enc->predictor[cl][k] = prd_save[cl][k];
            }

            i = 0;

            for (k = 0; k < enc->num_group; k++) {
                enc->th[cl][k] = th_save[cl][k];

                for (; i < enc->th[cl][k]; i++) {
                    enc->uquant[cl][i] = (char) k;
                }
            }
        }

        predict_region(enc, 0, 0, 0, 0, enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH]);
        cost = calc_cost(enc, 0, 0, 0, 0, enc->vu[HEIGHT], enc->vu[WIDTH], enc->ts[HEIGHT], enc->ts[WIDTH]);
        optimize_class(enc);
    }

    remove_emptyclass(enc);

    printf("2nd optimization --> Cost: %d (%d)", (int) cost, (int) side_cost);
    printf(" --> M = %d\n", enc->num_class);

    header = write_header(enc, prd_order, mi_prd_order, fp);

    if (forward_table != NULL) {
        header += encode_lookuptable(fp, forward_table, (int) pow(2, depth));
    }
    else {
        header += putbits(fp, 16, 0);
    }

    class_info = write_class(enc, fp);

    if (enc->f_huffman == 0) {
        enc->rc = rc_init();
    }

    class_info += encode_class(fp, enc);
    predictors = encode_predictor(fp, enc);
    thresholds = encode_threshold(fp, enc);
    errors = encode_image(fp, enc);

    free(class_save);
    free(prd_save);
    free(th_save);

    if (debug == 1) {
        debug_predictors(enc);
        debug_partition(enc, endianness);

        free(debug_path);
    }

    free_encoder(enc);

    safefree_lf4d(&lf);

    if (f_huffman == 1) {
        putbits(fp, 7, 0);    /* flush remaining bits */
    }

    fclose(fp);

    print_results(res, vu, mi_size, header, class_info, predictors, thresholds, errors);

    free(error);

    safefree((void **) &forward_table);

    elapse += cpu_time();

    fprintf(res, "\nCPU time: %.2f sec.\n\n", elapse);
    fclose(res);
    printf("\ncpu time: %.2f sec.\n", elapse);

    return (0);
}
