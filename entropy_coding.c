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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "mrp.h"

int *gen_hufflen(uint *hist, int size, int max_len) {
    int i, j, k, l, *len, *index, *bits, *link;

    len = (int *) alloc_mem(size * sizeof(int));
    index = (int *) alloc_mem(size * sizeof(int));
    bits = (int *) alloc_mem(size * sizeof(int));
    link = (int *) alloc_mem(size * sizeof(int));
    for (i = 0; i < size; i++) {
        len[i] = 0;
        index[i] = i;
        link[i] = -1;
    }
    /* sort in decreasing order of frequency */
    for (i = size - 1; i > 0; i--) {
        for (j = 0; j < i; j++) {
            if (hist[index[j]] < hist[index[j + 1]]) {
                k = index[j + 1];
                index[j + 1] = index[j];
                index[j] = k;
            }
        }
    }

    for (i = 0; i < size; i++) {
        bits[i] = index[i];    /* reserv a sorted index table */
    }

    for (i = size - 1; i > 0; i--) {
        k = index[i - 1];
        l = index[i];
        hist[k] += hist[l];
        len[k]++;

        while (link[k] >= 0) {
            k = link[k];
            len[k]++;
        }

        link[k] = l;
        len[l]++;

        while (link[l] >= 0) {
            l = link[l];
            len[l]++;
        }

        for (j = i - 1; j > 0; j--) {
            if (hist[index[j - 1]] < hist[index[j]]) {
                k = index[j];
                index[j] = index[j - 1];
                index[j - 1] = k;
            }
            else {
                break;
            }
        }
    }

    /* limit the maximum code length to max_len */
    for (i = 0; i < size; i++) {
        index[i] = bits[i];    /* restore the index table */
        bits[i] = 0;
    }

    for (i = 0; i < size; i++) {
        bits[len[i]]++;
    }

    for (i = size - 1; i > max_len; i--) {
        while (bits[i] > 0) {
            j = i - 2;

            while (bits[j] == 0) j--;

            bits[i] -= 2;
            bits[i - 1]++;
            bits[j + 1] += 2;
            bits[j]--;
        }
    }

    for (i = k = 0; i < size; i++) {
        for (j = 0; j < bits[i]; j++) {
            len[index[k++]] = i;
        }
    }
    free(link);
    free(bits);
    free(index);

    return (len);
}

void gen_huffcode(VLC *vlc) {
    int i, j, *idx, *len;
    uint k;

    vlc->index = idx = (int *) alloc_mem(vlc->size * sizeof(int));
    vlc->off = (int *) alloc_mem(vlc->max_len * sizeof(int));
    vlc->code = (uint *) alloc_mem(vlc->size * sizeof(int));
    len = vlc->len;

    /* sort in increasing order of code length */
    for (i = 0; i < vlc->size; i++) {
        idx[i] = i;
    }

    for (i = vlc->size - 1; i > 0; i--) {
        for (j = 0; j < i; j++) {
            if (len[idx[j]] > len[idx[j + 1]]) {
                k = (uint) (idx[j + 1]);
                idx[j + 1] = idx[j];
                idx[j] = k;
            }
        }
    }

    k = 0;

    for (j = 0; j < vlc->max_len; j++) {
        vlc->off[j] = -1;
    }

    j = len[idx[0]];

    for (i = 0; i < vlc->size; i++) {
        if (j < len[idx[i]]) {
            k <<= (len[idx[i]] - j);
            j = len[idx[i]];
        }

        vlc->code[idx[i]] = k++;
        vlc->off[j - 1] = i;
    }
}

VLC *make_vlc(uint *hist, int size, int max_len) {
    VLC *vlc;

    vlc = (VLC *) alloc_mem(sizeof(VLC));
    vlc->size = size;
    vlc->max_len = max_len;
    vlc->len = gen_hufflen(hist, size, max_len);
    gen_huffcode(vlc);
    return (vlc);
}

void free_vlc(VLC *vlc) {
    free(vlc->code);
    free(vlc->off);
    free(vlc->index);
    free(vlc->len);
    free(vlc);
}

VLC **init_vlcs(PMODEL ***pmodels, int num_group, int num_pmodel) {
    VLC **vlcs, *vlc;
    PMODEL *pm;
    int gr, k;

    vlcs = (VLC **) alloc_2d_array(num_group, num_pmodel, sizeof(VLC));
    for (gr = 0; gr < num_group; gr++) {
        for (k = 0; k < num_pmodel; k++) {
            vlc = &vlcs[gr][k];
            pm = pmodels[gr][k];
            vlc->size = pm->size;
            vlc->max_len = VLC_MAXLEN;
            vlc->len = gen_hufflen(pm->freq, pm->size, VLC_MAXLEN);
            gen_huffcode(vlc);
        }
    }
    return (vlcs);
}

/*
  Natural logarithm of the gamma function
  cf. "Numerical Recipes in C", 6.1
  http://www.ulib.org/webRoot/Books/Numerical_Recipes/bookcpdf.html
 */
double lngamma(double xx) {
    int j;
    double x, y, tmp, ser;
    double cof[6] = {
            76.18009172947146, -86.50532032941677,
            24.01409824083091, -1.231739572450155,
            0.1208650973866179e-2, -0.5395239384953e-5
    };

    y = x = xx;
    tmp = x + 5.5 - (x + 0.5) * log(x + 5.5);
    ser = 1.000000000190015;

    for (j = 0; j <= 5; j++) {
        ser += (cof[j] / ++y);
    }

    return (log(2.5066282746310005 * ser / x) - tmp);
}

void set_freqtable(PMODEL *pm, double *pdfsamp, int num_subpm, int num_pmodel, int center, int idx, double sigma) {
    double shape, beta, norm, sw, x;
    int i, j, n;

    if (center == 0) sigma *= 2.0;

    if (idx < 0) {
        shape = 2.0;
    }
    else {
        shape = 3.2 * (idx + 1) / (double) num_pmodel;
    }

    /* Generalized Gaussian distribution */
    beta = exp(0.5 * (lngamma(3.0 / shape) - lngamma(1.0 / shape))) / sigma;
    sw = 1.0 / (double) num_subpm;
    n = pm->size * num_subpm;
    center *= num_subpm;

    if (center == 0) {    /* one-sided distribution */
        for (i = 0; i < n; i++) {
            x = (double) i * sw;
            pdfsamp[i] = exp(-pow(beta * x, shape));
        }
    }
    else {
        for (i = center; i < n; i++) {
            x = (i - (double) center + 0.5) * sw;
            pdfsamp[i + 1] = exp(-pow(beta * x, shape));
        }

        for (i = 0; i <= center; i++) {
            pdfsamp[center - i] = pdfsamp[center + i + 1];
        }

        for (i = 0; i < n; i++) {
            if (i == center) {
                pdfsamp[i] = (2.0 + pdfsamp[i] + pdfsamp[i + 1]) / 2.0;
            }
            else {
                pdfsamp[i] = pdfsamp[i] + pdfsamp[i + 1];
            }
        }
    }

    for (j = 0; j < num_subpm; j++) {
        norm = 0.0;
        for (i = 0; i < pm->size; i++) {
            norm += pdfsamp[i * num_subpm + j];
        }

        norm = (double) (MAX_TOTFREQ - pm->size * MIN_FREQ) / norm;
        norm += 1E-8;    /* to avoid machine dependent rounding errors */
        pm->cumfreq[0] = 0;

        for (i = 0; i < pm->size; i++) {
            pm->freq[i] = (uint) (norm * pdfsamp[i * num_subpm + j] + MIN_FREQ);
            pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
        }

        pm++;
    }
}

/*---------------------------- init_pmodels ------------------------*
 |  Function init_pmodels
 |
 |  Purpose:  Initializes the probability models
 |
 |  Parameters:
 |		num_group		--> Number of groups ???? (IN)
 |		num_pmodel		--> Number of probability models (IN)
 |		pm_accuracy		--> Probability model accuracy (IN)
 |		pm_idx			--> ?? --> NULL (IN)
 |		sigma			--> Predefined vectors (IN)
 |		size			--> Image maximum value + 1 (IN)
 |
 |  Returns:  PMODEL*** --> returns a PMODEL type structure
 *-------------------------------------------------------------------*/
PMODEL ***init_pmodels(int num_group, int num_pmodel, int pm_accuracy, int *pm_idx, double *sigma, int size) {
    PMODEL ***pmodels, *pmbuf, *pm;
    int gr, i, j, num_subpm, num_pm, ssize, idx;
    double *pdfsamp;

    //??
    if (pm_accuracy < 0) {
        num_subpm = 1;
        ssize = 1;
    }
    else {
        num_subpm = 1 << pm_accuracy;
        ssize = size;
        size = size + ssize - 1;
    }

    //Defines the number of probability models
    num_pm = (pm_idx != NULL) ? 1 : num_pmodel;

    pmodels = (PMODEL ***) alloc_2d_array(num_group, num_pm, sizeof(PMODEL *));
    pmbuf = (PMODEL *) alloc_mem(num_group * num_pm * num_subpm * sizeof(PMODEL));

    //Vector initializations
    for (gr = 0; gr < num_group; gr++) {
        for (i = 0; i < num_pm; i++) {
            pmodels[gr][i] = pmbuf;

            for (j = 0; j < num_subpm; j++) {
                pm = pmbuf++;
                pm->id = i;
                pm->size = size;
                pm->freq = (uint *) alloc_mem((size * 2 + 1) * sizeof(uint));
                pm->cumfreq = &pm->freq[size];

                if (pm_idx == NULL) {
                    pm->cost = (float *) alloc_mem((size + ssize) * sizeof(float));
                    pm->subcost = &pm->cost[size];
                }
            }
        }
    }

    pdfsamp = alloc_mem((size * num_subpm + 1) * sizeof(double));

    for (gr = 0; gr < num_group; gr++) {
        for (i = 0; i < num_pm; i++) {
            if (pm_idx != NULL) {
                idx = pm_idx[gr];
            }
            else if (num_pm > 1) {
                idx = i;
            }
            else {
                idx = -1;
            }

            set_freqtable(pmodels[gr][i], pdfsamp, num_subpm, num_pmodel, ssize - 1, idx, sigma[gr]);
        }
    }

    free(pdfsamp);

    return (pmodels);
}

/* probability model for coefficients and thresholds */
void set_spmodel(PMODEL *pm, int size, int m) {
    int i, sum;
    double p;

    pm->size = size;

    if (m >= 0) {
        p = 1.0 / (double) (1 << (m % 8));
        sum = 0;

        for (i = 0; i < pm->size; i++) {
            pm->freq[i] = (uint) (exp(-p * i) * (1 << 10));

            if (pm->freq[i] == 0) pm->freq[i]++;

            sum += pm->freq[i];
        }

        if (m & 8) pm->freq[0] = (sum - pm->freq[0]);    /* weight for zero */
    }
    else {
        for (i = 0; i < pm->size; i++) {
            pm->freq[i] = 1;
        }
    }

    pm->cumfreq[0] = 0;

    for (i = 0; i < pm->size; i++) {
        pm->cumfreq[i + 1] = pm->cumfreq[i] + pm->freq[i];
    }
}