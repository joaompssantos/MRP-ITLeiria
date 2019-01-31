#ifndef MRP_H
#define MRP_H

#define HAVE_CLOCK
#define HAVE_64BIT_INTEGER
#define MAGIC_NUMBER        ('M' << 8) + 'R'
#define BANNER              "\nIT - Leiria: Minimum Rate Predictors\nencmrp/decmrp version %s (%s)"
#define uint                unsigned int
#define img_t               unsigned short
#define cost_t              double

// Entropy coding
#ifdef HAVE_64BIT_INTEGER
#  define RANGE_SIZE 64
#  if defined(_MSC_VER) || defined(__BORLANDC__)
#    define range_t unsigned __int64
#  else
#    define range_t unsigned long long
#  endif
#  define MAX_TOTFREQ (1 << 20)	/* must be < RANGE_BOT */
#else
#  define RANGE_SIZE 32
#  define range_t unsigned int
#  define MAX_TOTFREQ (1 << 14)	/* must be < RANGE_BOT */
#endif
#define RANGE_TOP  ((range_t)1 << (RANGE_SIZE - 8))
#define RANGE_BOT  ((range_t)1 << (RANGE_SIZE - 16))

#define NUM_PMODEL          16
#define PM_ACCURACY         3
#define MIN_FREQ            1
#define MAX_SYMBOL          1024        /* must be >> MAX_UPARA */
#define PMCLASS_MAX         16
#define PMCLASS_LEVEL       32
#define VLC_MAXLEN          32

// MRP Common
#define FRAMES              1
#define DEPTH               8
#define QUADTREE_DEPTH      4
#define BASE_BSIZE          8
#define MAX_BSIZE           32
#define MIN_BSIZE           (MAX_BSIZE >> QUADTREE_DEPTH)
#define MAX_CLASS           63
#define NUM_CLASS           -1
#define NUM_GROUP           16
#define INTRA_PRD_ORDER     -1
#define INTER_PRD_ORDER     1
#define COEF_PRECISION      6
#define MAX_COEF            (2 << COEF_PRECISION)
#define MAX_UPARA           512
#define MAX_ITERATION       100
#define EXTRA_ITERATION     10
#define OPT_SIDEINFO

// Endianness
#define LITTLE_ENDIANNESS   0
#define BIG_ENDIANNESS      1

// Chroma format
#define GRAY                0
#define S444                1
#define S422                2
#define S411                3
#define S420                4

// Image
typedef struct {
    int height;
    int width;
    int maxval;
    img_t **val;
} IMAGE;

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

// Point
typedef struct {
    int y, x;
} POINT;

// Encoder
typedef struct {
    int height; // Image height
    int width; // Image width
    int delta; // Distance between the reference frame and the current one
    int depth; // Bit depth of the input image/sequence
    int maxval; // Image maximum value
    int num_class; // Number of classes to use (number of different predictors)
    int num_group; // Number of pixels that are taken into account to form the context for entropy coding of the residue (= 16)
    int prd_order; // Order of the predictors (number of pixels to use)
    int back_prd_order; // Order of the predictors (number of pixels to use) in the previous frame
    int for_prd_order; // Order of the predictors (number of pixels to use) in the next frame
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
    int **encval; // Pointer to the error (in set_cost_model()) and later to the image real values (in set_cost_rate()).
    int ***err; // Matrix that keeps residue values after the prediction.
    int ***org; // Original video/image
    int *ctx_weight; // Keeps the weights used for the residue encoding context.
    int ***roff; // Auxiliary structure used to scan the image for neighboring pixels.
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
} ENCODER;

// Decoder
typedef struct {
    int version;
    int height;
    int width;
    int maxval;
    int frames;
    int delta;
    int depth; // Bit depth of the input image/sequence
    int num_comp;
    int num_class;
    int num_group;
    int prd_order;
    int back_prd_order; // Order of the predictors (number of pixels to use) in the previous frame
    int for_prd_order; // Order of the predictors (number of pixels to use) in the next frame
    int num_pmodel;
    int pm_accuracy;
    int maxprd;
    int coef_precision;
    int f_huffman;
    int quadtree_depth;
    int **predictor;
    int ***err;
    int *ctx_weight;
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
} DECODER;

/* common.c */
void safefree(void **);
void safefree_yuv(IMAGE **);
FILE *fileopen(char *, char *);
void *alloc_mem(size_t);
void **alloc_2d_array(int, int, int);
void ***alloc_3d_array(int, int, int, int);
IMAGE *alloc_image(int, int, int);
IMAGE *copy_yuv(IMAGE*);
unsigned short reverse_endianness (unsigned short, int);
void write_yuv(IMAGE*, char*, int, int);
IMAGE *read_yuv(char *, int, int, int, int, int, int);
int *gen_hufflen(uint *, int, int);
void gen_huffcode(VLC *);
VLC *make_vlc(uint *, int, int);
VLC ***init_vlc(int, int, int *, int);
void free_vlc(VLC *);
VLC **init_vlcs(PMODEL ***, int, int);
PMODEL ***init_pmodels(int, int, int, int *, double *, int);
void set_spmodel(PMODEL *, int, int);
int *init_ctx_weight(int, int, int, int);
int e2E(int, int, int, int);
int E2e(int, int, int, int);
void mtf_classlabel(char **, int *, int, int, int, int, int);
double cpu_time(void);
int **select_bref(int);
char* cat_str(char*, char*, int);
char *int2bin(int, int);

/* rc.c */
RANGECODER *rc_init(void);
void rc_encode(FILE *, RANGECODER *, uint, uint, uint);
void rc_finishenc(FILE *, RANGECODER *);
int rc_decode(FILE *, RANGECODER *, PMODEL *, int, int);
void rc_startdec(FILE *, RANGECODER *);

#endif //MRP_H