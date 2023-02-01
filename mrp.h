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

#ifndef MRP_H
#define MRP_H

#define HAVE_CLOCK
#define HAVE_64BIT_INTEGER
#define MAGIC_NUMBER        ('M' << 8) + 'R'
#define BANNER              "\nIT - Leiria: Hierarchical Minimum Rate Predictors\nencmrp/decmrp version %s (%s)"
#define uint                unsigned int
#define img_t               unsigned short
#define cost_t              double

// Entropy coding
#ifdef HAVE_64BIT_INTEGER
#  define RANGE_SIZE        64
#  if defined(_MSC_VER) || defined(__BORLANDC__)
#    define range_t         unsigned __int64
#  else
#    define range_t         unsigned long long
#  endif
#  define MAX_TOTFREQ       (1 << 20)    /* must be < RANGE_BOT */
#else
#  define RANGE_SIZE        32
#  define range_t           unsigned int
#  define MAX_TOTFREQ       (1 << 14)	/* must be < RANGE_BOT */
#endif
#define RANGE_TOP           ((range_t) 1 << (RANGE_SIZE - 8))
#define RANGE_BOT           ((range_t) 1 << (RANGE_SIZE - 16))

#define NUM_PMODEL          16
#define PM_ACCURACY         3
#define MIN_FREQ            1
#define MAX_SYMBOL          1024        /* must be >> MAX_UPARA */
#define PMCLASS_MAX         16
#define PMCLASS_LEVEL       32
#define VLC_MAXLEN          32

// MRP Common
#define DEPTH               8
#define QUADTREE_DEPTH      4
#define BASE_BSIZE          8
#define MAX_BSIZE           32
#define MIN_BSIZE           (MAX_BSIZE >> QUADTREE_DEPTH)
#define MAX_CLASS           63
#define NUM_CLASS           -1
#define NUM_GROUP           16
#define INTRA_PRD_ORDER     -1
#define SAI_PRD_ORDER        1
#define COEF_PRECISION      6
#define MAX_COEF            (2 << COEF_PRECISION)
#define MAX_UPARA           512
#define MAX_ITERATION       100
#define EXTRA_ITERATION     10
#define OPT_SIDEINFO

#define HEIGHT              0
#define WIDTH               1

// Endianness
#define LITTLE_ENDIANNESS   0
#define BIG_ENDIANNESS      1

// Chroma format
#define GRAY                0
#define S444                1
#define S422                2
#define S411                3
#define S420                4

// Data format
#define MIA                 0
#define PVS                 1
#define SAI                 2

// Hierarchical coding
#define SAI_REFERENCES             4
#define SAI_DISTANCE_THRESHOLD     8

// Disparity estimation
#define DISP_BSIZE                 8
#define MAX_DISPARITY              8

// Entropy coding
typedef struct {
    int size;
    int id;
    uint *freq;
    uint *cumfreq;
    float *cost;
    float *subcost;
} PMODEL;

typedef struct {
    int size;
    int max_len;
    int *len;
    int *index;
    int *off;
    uint *code;
} VLC;

typedef struct {
    range_t low;
    range_t code;
    range_t range;
} RANGECODER;

// 4DLF
typedef struct {
    int t; // View height: t
    int s; // View width: s
    int v; // Views array height: v
    int u; // Views array width: u
    int maxval;
    img_t ****val;
} LF4D;

// Point
typedef struct {
    int y, x;
} POINT;

// Store the SAI number and its distance
typedef struct {
    int sai;
    int coordinates[2];
    double distance;
} SAIDISTANCE;

// Encoder
typedef struct {
    int ts[2]; // View dimensions

    int no_refsai; // Number of reference SAIs available
    int delta; // Distance between the reference frame and the current one
    int depth; // Bit depth of the input image/sequence
    int maxval; // Image maximum value
    int num_class; // Number of classes to use (number of different predictors)
    int num_group; // Number of pixels that are taken into account to form the context for entropy coding of the residue (= 16)

    int full_prd_order;
    int prd_order; // Order of the predictors (number of pixels to use)
    int sai_prd_order; // Order of the predictors (number of pixels to use) in neighbouring micro images

    int coef_precision; // Precision of the coefficients
    int num_pmodel; // Number of probability models
    int pm_accuracy; // Probability model accuracy
    int maxprd; // Maximum prediction value allowed, depends on the coef_precision
    int f_huffman; // Huffman coding flag
    int quadtree_depth; // Quadtree depth on/off
    int optimize_loop; // First or second optimization loop
    int **predictor; // Stores the prediction coefficients for all classes
    int **th; // Indicates the threshold values used to quantize the weighted sum of neighboring residues, to be used on the context coding.
    int **upara; // Context for the residue encoding, without the quantization.
    int **prd; // Prediction value for each pixel
    int ***err; // Matrix that keeps residue values after the prediction.
    int ***org; // Original video/image
    int **encval; // Pointer to the error (in set_cost_model()) and later to the image real values (in set_cost_rate()).

    int *intra_ctx_weight; // Keeps the weights used for the residue encoding context.
    int *sai_ref_ctx_weight[SAI_REFERENCES]; // Keeps the weights used for the residue encoding context.

    int ***intra_roff; // Auxiliary structure used to scan the image for neighboring pixels.
    int ***sai_ref_roff[SAI_REFERENCES]; // Auxiliary structure used to scan the image for neighboring pixels.
    int ***disparity_vectors[SAI_REFERENCES]; // Saves the disparity vectors for each reference image.

    int qtctx[QUADTREE_DEPTH << 3]; // Frequency counter of the segmentation flag context for the whole image.
    char **qtmap[QUADTREE_DEPTH]; // Segmentation flags for the quadtree partitioning of the image's prediction.

    char **class; // Keeps the class of each pixel
    char **group; // keeps the group, i.e. the quantized context used for residue entropy coding.

    char **uquant; // Table used for the quantification of the context variable u.
    int etype; // Error type for the error conversion table
    int **econv; // Table used for the conversion of the error.
    img_t *bconv; // Structure used to convert the prediction to a pointer which indicates the position in the probability vector structure of the prediction error.
    img_t *fconv; // Structure used to fine tune the probability value, given the probability model accuracy.
    PMODEL ***pmodels; // Structure that holds all probability models.
    PMODEL **pmlist; // List of pointer for the probability model for each group.
    PMODEL spm; // Probability model structure used for the side information (classes, thresholds, coefficients).
    VLC **vlcs; // Structure for the Huffman encoder.
    RANGECODER *rc; // Structure for the Range encode.
    double *sigma; // A vector that defines the variance that will be used for the generalized gaussian, for each context.
    int *mtfbuf; // Buffer used for context determination for class encoding.
    int *coef_m; // Structure that indicates the context used for arithmetic encoding of the coefficients of the prediction filters.
    cost_t **coef_cost; // Structure used to keep the cost of the coefficients.
    cost_t *th_cost; // Structure used to keep the cost of the thresholds.
    cost_t *class_cost; // Array with the cost of each class.
    cost_t qtflag_cost[QUADTREE_DEPTH << 3]; // Structure for the cost of the segmentation flags.

    char *debug_path;
} ENCODER;

// Decoder
typedef struct {
    int ts[2]; // View dimensions

    int no_refsai; // Number of reference SAIs available
    int maxval;
    int delta;
    int depth; // Bit depth of the input image/sequence
    int num_comp;
    int num_class;
    int num_group;

    int full_prd_order;
    int prd_order;
    int sai_prd_order; // Order of the predictors (number of pixels to use) in neighbouring micro images

    int ***intra_roff; // Auxiliary structure used to scan the image for neighboring pixels.
    int ***sai_ref_roff[SAI_REFERENCES]; // Auxiliary structure used to scan the image for neighboring pixels.
    int ***disparity_vectors[SAI_REFERENCES]; // Saves the disparity vectors for each reference image.
    
    int num_pmodel;
    int pm_accuracy;
    int maxprd;
    int coef_precision;
    int f_huffman;
    int quadtree_depth;
    int **predictor;
    int ***err;
    img_t ***org; // Original video/image

    int *intra_ctx_weight; // Keeps the weights used for the residue encoding context.
    int *sai_ref_ctx_weight[SAI_REFERENCES]; // Keeps the weights used for the residue encoding context.

    char **qtmap[QUADTREE_DEPTH];
    char **class;

    char **uquant;
    int *pm_idx;
    PMODEL ***pmodels;
    PMODEL spm;
    VLC **vlcs;
    RANGECODER *rc;
    double *sigma;
    int *mtfbuf;

    char *debug_path;
} DECODER;

/* common.c */
void safefree(void **);

void safefree_lf4d(LF4D **);

FILE *fileopen(char *, char *);

void *alloc_mem(size_t);

void **alloc_2d_array(size_t, size_t, size_t);

void ***alloc_3d_array(size_t, size_t, size_t, size_t);

void ****alloc_4d_array(size_t, size_t, size_t, size_t, size_t);

LF4D *alloc_lf4d(int, int, int, int, int);

LF4D *copy_yuv(LF4D *);

unsigned short reverse_endianness(unsigned short, int);

void write_yuv(LF4D *, char *, int, int, int);

LF4D *read_yuv(char *, int, int, int, int, int, int, int);

int ***init_ref_offset(int[2], int, int, int ***, int);

void free_ref_offset(const int[2], int, int, int ***);

int *init_ctx_weight(int, int, int);

int e2E(int, int, int, int);

int E2e(int, int, int, int);

void mtf_classlabel(char **, int *, int, int, int, int, int);

double cpu_time(void);

char *cat_str(char *, char *, int);

char *int2bin(int, int);

void frame2coordinates(int *, int, int);

int compare_distance(const void *, const void *);

double median(int [], int);

int get_random_access_region(int, const int [2], int, int, int [2]);


/* Huffman */
int *gen_hufflen(uint *, int, int);

void gen_huffcode(VLC *);

VLC *make_vlc(uint *, int, int);

VLC ***init_vlc(int, int, int *, int);

void free_vlc(VLC *);

VLC **init_vlcs(PMODEL ***, int, int);

/* Probability models */
PMODEL ***init_pmodels(int, int, int, int *, double *, int);

void set_spmodel(PMODEL *, int, int);

void update_pmodel(PMODEL *, int, int);

/* rc.c */
RANGECODER *rc_init(void);

void rc_encode(FILE *, RANGECODER *, uint, uint, uint);

void rc_finishenc(FILE *, RANGECODER *);

int rc_decode(FILE *, RANGECODER *, PMODEL *, int, int);

void rc_startdec(FILE *, RANGECODER *);

#endif //MRP_H
