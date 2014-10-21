/*=============================================================================
 |   Assignment:  ASSIGNMENT NUMBER AND TITLE
 |
 |       Author:  STUDENT'S NAME HERE
 |     Language:  NAME OF LANGUAGE IN WHICH THE PROGRAM IS WRITTEN AND THE
 |                      NAME OF THE COMPILER USED TO COMPILE IT WHEN IT
 |                      WAS TESTED
 |   To Compile:  EXPLAIN HOW TO COMPILE THIS PROGRAM
 |
 |        Class:  NAME AND TITLE OF THE CLASS FOR WHICH THIS PROGRAM WAS
 |                      WRITTEN
 |   Instructor:  NAME OF YOUR COURSE'S INSTRUCTOR
 |     Due Date:  DATE AND TIME THAT THIS PROGRAM IS/WAS DUE TO BE SUBMITTED
 |
 +-----------------------------------------------------------------------------
 |
 |  Description:  DESCRIBE THE PROBLEM THAT THIS PROGRAM WAS WRITTEN TO
 |      SOLVE.
 |
 |        Input:  DESCRIBE THE INPUT THAT THE PROGRAM REQUIRES.
 |
 |       Output:  DESCRIBE THE OUTPUT THAT THE PROGRAM PRODUCES.
 |
 |    Algorithm:  OUTLINE THE APPROACH USED BY THE PROGRAM TO SOLVE THE
 |      PROBLEM.
 |
 |   Required Features Not Included:  DESCRIBE HERE ANY REQUIREMENTS OF
 |      THE ASSIGNMENT THAT THE PROGRAM DOES NOT ATTEMPT TO SOLVE.
 |
 |   Known Bugs:  IF THE PROGRAM DOES NOT FUNCTION CORRECTLY IN SOME
 |      SITUATIONS, DESCRIBE THE SITUATIONS AND PROBLEMS HERE.
 |
 *===========================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include "mrp.h"
#include <omp.h>

extern POINT dyx[];
extern POINT idyx[];
extern double sigma_h[], sigma_a[];
extern double qtree_prob[];

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

int pixel_wise_predictor(ENCODER *enc){
	int f, i, j, aux;
	int resultado = 0;
	IMAGE *img = alloc_image(enc->width, enc->height, enc->frames, 255);
	IMAGE *img2 = alloc_image(enc->width, enc->height, enc->frames, 255);

	for(f = 0; f < enc->frames; f++){
		for(i = 0; i < enc->height; i++){
			for(j = 0; j < enc->width; j++){
				if(f == 0){
					img->val[f][i][j] = enc->org[f][i][j];
					img2->val[f][i][j] = enc->org[f][i][j];
				}
				else{
					aux = enc->org[f][i][j] - enc->org[f - 1][i][j];
					if(aux > 128){
						aux = 128;
						resultado = 1;
					}
					else if(aux < -127){
						aux = -127;
						resultado = 1;
					}
					img->val[f][i][j] =  aux + 127;
					img2->val[f][i][j] = enc->org[f - 1][i][j];
				}
			}
		}
	}

	//write_yuv(img2, "pixel_wise_predictor.yuv");
	
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

int pixel_wise_search_predictor(ENCODER *enc, int window){
	int f, i, j, k, z, aux, dif = 512, dif_aux;
	int pos_x, pos_y;
	int resultado = 0;
	IMAGE *img = alloc_image(enc->width, enc->height, enc->frames, 255);
	IMAGE *img2 = alloc_image(enc->width, enc->height, enc->frames, 255);

	for(f = 0; f < enc->frames; f++){
		for(i = 0; i < enc->height; i++){
			for(j = 0; j < enc->width; j++){
				if(f == 0){
					img->val[f][i][j] = enc->org[f][i][j];
					img2->val[f][i][j] = enc->org[f][i][j];
				}
				else if(i < window + 1 || j < window + 1 || i > enc->height - window - 1|| j > enc->width - window - 1){
					img->val[f][i][j] = enc->org[f][i][j];
					img2->val[f][i][j] = enc->org[f][i][j];
				}
				else{
					dif = 512;
// 					pos_y = i;
// 					pos_x = j;
// 				  
// 					dif = enc->org[f][i][j - 0] - enc->org[f - 1][i][j - 0];
// 				  
// 					if(dif != 0){
						for(k = i - window; k <= i + window; k++){
							for(z = j - window; z <= j + window; z++){
								dif_aux = enc->org[f][i][j - 0] - enc->org[f - 1][k][z - 0];
								
								if(abs(dif_aux) < abs(dif)){
									//if(dif_aux != 0){printf("teste: %d %d\n", dif_aux, dif);getchar();}
									pos_y = k;
									pos_x = z;
									
									dif = dif_aux;
								}
							}
						}
// 					}
					
//  					if(dif != 0) {
// 						printf("Dif: %d, org: %d, novo: %d, posI: %d %d, posF: %d %d\n", dif, enc->org[f][i][j - 0], enc->org[f - 1][pos_y][pos_x], i, j, pos_y, pos_x);
// 
// 						printf("Valor a prever = %d\nMatriz inicial\n\n", enc->org[f][i][j]);
// 						for(k = i - window; k <= i + window; k++){
// 							for(z = j - window; z <= j + window; z++){
// 								printf("%3d ", enc->org[f][k][z]);
// 							}
// 							printf("\n");
// 						}
// 						printf("Matriz Outra:\n\n");
// 						for(k = i - window; k <= i + window; k++){
// 							for(z = j - window; z <= j + window; z++){
// 								printf("%3d ", enc->org[f - 1][k][z]);
// 							}
// 							printf("\n");
// 						}
// 						getchar();
// 					}
					//printf("diferenca = %d\n", dif);
					aux = enc->org[f][i][j] - enc->org[f - 1][pos_y][pos_x];
					if(aux > 128){
						aux = 128;
						resultado = 1;
					}
					else if(aux < -127){
						aux = -127;
						resultado = 1;
					}
					
					img->val[f][i][j] =  aux + 127;
					img2->val[f][i][j] = enc->org[f - 1][pos_y][pos_x];
				}
			}
		}
	}

	write_yuv(img2, "pixel_wise_search_predictor.yuv");
	
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

int gap_predictor(ENCODER *enc){
	int f, i, j, aux, frame, resultado = 0;
	int dh, dv, d45, d135;
	IMAGE *img = alloc_image(enc->width, enc->height, enc->frames, 255);

	for(f = 0; f < enc->frames; f++){
		for(i = 0; i < enc->height; i++){
			for(j = 0; j < enc->width; j++){
				if(f == 0){
					img->val[f][i][j] = enc->org[f][i][j];
				}
				else if(i < 2 || j < 2 || i > enc->height - 2 || j > enc->width - 2){
					img->val[f][i][j] = enc->org[f][i][j];
				}
				else{
					frame = f - 1;

					dh = abs(enc->org[frame][i][j - 1] - enc->org[frame][i][j - 2]) + abs(enc->org[frame][i - 1][j] - enc->org[frame][i - 1][j - 1]) + abs(enc->org[frame][i - 1][j + 1] - enc->org[frame][i - 1][j]);
					dv = abs(enc->org[frame][i][j - 1] - enc->org[frame][i - 1][j - 1]) + abs(enc->org[frame][i - 1][j] - enc->org[frame][i - 2][j]) + abs(enc->org[frame][i - 1][j + 1] - enc->org[frame][i - 2][j + 1]);
					d45 = abs(enc->org[frame][i][j - 1] - enc->org[frame][i - 1][j - 2]) + abs(enc->org[frame][i - 1][j - 1] - enc->org[frame][i - 2][j - 2]) + abs(enc->org[frame][i - 1][j] - enc->org[frame][i - 2][j - 1]);
					d135 = abs(enc->org[frame][i - 1][j + 1] - enc->org[frame][i - 2][j + 2]) + abs(enc->org[frame][i - 1][j] - enc->org[frame][i - 2][j + 1]) + abs(enc->org[frame][i - 1][j] - enc->org[frame][i][j - 1]);

					if(dv + dh > 32){ //Sharp Edge
						aux = (dv * enc->org[frame][i][j - 1] + dh * enc->org[frame][i - 1][j]) / (dv + dh) + (enc->org[frame][i - 1][j + 1] - enc->org[frame][i - 1][j - 1]) / 8;
					}
					else if(dv -dh > 12){ //Horizontal Edge
						aux = (2 * enc->org[frame][i][j - 1] + enc->org[frame][i - 1][j]) / 3 + (enc->org[frame][i - 1][j + 1] - enc->org[frame][i - 1][j - 1]) / 8;
					}
					else if(dh - dv > 12){ //Vertical Edge
						aux = (enc->org[frame][i][j - 1] + 2 * enc->org[frame][i - 1][j]) / 3 + (enc->org[frame][i - 1][j + 1] - enc->org[frame][i - 1][j - 1]) / 8;
					}
					else{ //Smooth Area
						aux = (enc->org[frame][i][j - 1] + enc->org[frame][i - 1][j]) / 2 + (enc->org[frame][i - 1][j + 1] - enc->org[frame][i - 1][j - 1]) / 8;
					}

					if(d45 + d135 > 32){ //Sharp 135-deg diagonal-edge
						aux = aux + (enc->org[frame][i - 1][j + 1] - enc->org[frame][i - 1][j - 1]) / 8;
					}
					else if(d45 + d135 < 16){ //135-deg diagonal-edge
						aux = aux + (enc->org[frame][i - 1][j + 1] - enc->org[frame][i - 1][j - 1]) / 16;
					}
					else if(d135 - d45 > 32){ //Sharp 45-deg diagonal edge
						aux = aux + (enc->org[frame][i - 1][j - 1] - enc->org[frame][i - 1][j + 1]) / 8;
					}
					else if(d135 - d45 > 16){ //45-deg diagonal edge
						aux = aux + (enc->org[frame][i - 1][j - 1] - enc->org[frame][i - 1][j + 1]) / 16;
					}

					aux = aux - enc->org[f][i][j];
					if(aux > 128){
						aux = 128;
						resultado = 1;
					}
					else if(aux < -127){
						aux = -127;
						resultado = 1;
					}
					img->val[f][i][j] =  aux + 127;
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

	write_yuv(img, "gap_pred.yuv");
	free(img);
	return resultado;
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
 |		frames		--> Frames of the video (IN)
 |
 |  Returns:  IMAGE* --> returns a video type structure
 *-------------------------------------------------------------------*/
IMAGE *read_yuv(char *filename, int height, int width, int frames){
	int f, i, j;
	IMAGE *img;
	FILE *fp;

	//Open file
	fp = fileopen(filename, "rb");

	// Check if image dimensions are correct (It has to be multiple of BASE_BSIZE)
	if ((width % BASE_BSIZE) || (height % BASE_BSIZE)) {
		fprintf(stderr, "Image width and height must be multiples of %d!\n", BASE_BSIZE);
		exit(1);
	}

	// Image allocation
	img = alloc_image(width, height, frames, 255);
	for(f = 0; f < img->frames; f++){
		for (i = 0; i < img->height; i++) {
			for (j = 0; j < img->width; j++) {
				img->val[f][i][j] = (img_t)fgetc(fp);
			}
		}

		for (i = 0; i < img->height / 2; i++) {
			for (j = 0; j < img->width / 2; j++) {
				fgetc(fp);
			}
		}

		for (i = 0; i < img->height / 2; i++) {
			for (j = 0; j < img->width / 2; j++) {
				fgetc(fp);
			}
		}
	}

	fclose(fp);
	return (img);
}

/*-----------------------_--- init_ref_offset -----------------------*
 |  Function init_ref_offset
 |
 |  Purpose:  Calculates the offset of a reference to a given pixel
 |
 |  Parameters:
 |		img				--> Image structure (IN)
 |		prd_order		--> Number of intra reference pixels (IN)
 |		inter_prd_order	--> Number of inter reference pixels (IN)
 |
 |  Returns:  int*** --> Array with the references offset
 *-------------------------------------------------------------------*/
int ***init_ref_offset(IMAGE *img, int prd_order, int inter_prd_order){
	int ***roff, *ptr;
	int x, y, dx, dy, k;
	int order, min_dx, max_dx, min_dy;
	int imin_dx, imax_dx, imin_dy, imax_dy;
	int min_abs_dx, max_abs_dx, min_abs_dy, max_abs_dy;

	min_dx = max_dx = min_dy = 0;
	imin_dx = imax_dx = imin_dy = imax_dy = 0;
	order = (prd_order > NUM_UPELS)? prd_order : NUM_UPELS;

	//Values to check for special cases
	for (k = 0; k < order; k++) {
		dy = dyx[k].y;
		dx = dyx[k].x;
		if (dy < min_dy) min_dy = dy;
		if (dx < min_dx) min_dx = dx;
		if (dx > max_dx) max_dx = dx;
	}
	for (k = 0; k < inter_prd_order; k++) {
		dy = idyx[k].y;
		dx = idyx[k].x;

		if (dy < imin_dy) imin_dy = dy;
		if (dy > imax_dy) imax_dy = dy;
		if (dx < imin_dx) imin_dx = dx;
		if (dx > imax_dx) imax_dx = dx;
	}

	max_abs_dy = imax_dy;
	if(min_dy < imin_dy){
		min_abs_dy = min_dy;
	}
	else{
		min_abs_dy = imin_dy;
	}


	if(max_dx < imax_dx){
		max_abs_dx = imax_dx;
	}
	else{
		max_abs_dx = max_dx;
	}
	if(min_dx < imin_dx){
		min_abs_dx = min_dx;
	}
	else{
		min_abs_dx = imin_dx;
	}

	roff = (int ***)alloc_2d_array(img->height, img->width, sizeof(int *));
	//ptr = (int *)alloc_mem((1 - min_abs_dy + max_abs_dy) * (1 + max_abs_dx - min_abs_dx) * (order + inter_prd_order) * sizeof(int));

	//Cycle that runs for all the pixels
	for (y = 0; y < img->height; y++) {
		for (x = 0; x < img->width; x++) {
			ptr = (int *) alloc_mem((order + inter_prd_order) * sizeof(int));
			//Conditions to check which references are available for each pixel
			if (y == 0) {
				if (x == 0) {
					roff[y][x] = ptr;
					dx = 0;
					dy = img->height;

					for (k = 0; k < order; k++) {
						*ptr++ = dy * img->width + dx; //Points to a line filled with 128 (if max_val = 256)
					}
				}
				else if (x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
					roff[y][x] = ptr;
					dy = 0;

					for (k = 0; k < order; k++) {
						dx = dyx[k].x;

						if (x + dx < 0) dx = -x;
						else if (dx >= 0) dx = -1;

						*ptr++ = dy * img->width + dx;
					}
				}
				else {
					roff[y][x] = roff[y][x - 1];
				}
			}
			else if (y + min_abs_dy <= 0) {
				if (x == 0) {
					roff[y][x] = ptr;

					for (k = 0; k < order; k++) {
						dy = dyx[k].y;

						if (y + dy < 0) dy = -y;
						else if (dy >= 0) dy = -1;

						dx = dyx[k].x;

						if (x + dx < 0) dx = -x;

						*ptr++ = dy * img->width + dx;
					}
				}
				else if (x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
					roff[y][x] = ptr;

					for (k = 0; k < order; k++) {
						dy = dyx[k].y;

						if (y + dy < 0) dy = -y;
						dx = dyx[k].x;

						if (x + dx < 0) dx = -x;
						else if (x + dx >= img->width) {
							dx = img->width - x - 1;
						}

						*ptr++ = dy * img->width + dx;
					}
				}
				else {
					roff[y][x] = roff[y][x - 1];
				}
			}
			else if (y + max_abs_dy >= img->height) {
				roff[y][x] = ptr;

				for (k = 0; k < order; k++) {
					*ptr++ = roff[y - 1][x][k];
				}
			}
			else {
				roff[y][x] = roff[y - 1][x];
			}

			//Inter reference offset
			if(inter_prd_order != 0){
				int base = -img->width * (img->height + 1);

				if (y == 0) {
					if (x == 0) {
						*ptr++ = base;

						for (k = 0; k < inter_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if(y + dy < 0 || x + dx < 0){
								*ptr++ = base;
							}
							else{
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else if (x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
						*ptr++ = base;

						for (k = 0; k < inter_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if(y + dy < 0 || x + dx < 0 || x + dx >= img->width){
								*ptr++ = base;
							}
							else{
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else {
						roff[y][x] = roff[y][x - 1];
					}
				}
				else if (y + min_abs_dy <= 0) {
					if (x == 0) {
						*ptr++ = base;

						for (k = 0; k < inter_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if(y + dy < 0 || x + dx < 0){
								*ptr++ = base;
							}
							else{
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else if (x + min_abs_dx <= 0 || x + max_abs_dx >= img->width) {
						*ptr++ = base;

						for (k = 0; k < inter_prd_order - 1; k++) {
							dy = idyx[k].y;
							dx = idyx[k].x;

							if(x + dx < 0 || x + dx >= img->width){
								*ptr++ = base;
							}
							else{
								*ptr++ = dy * img->width + dx + base;
							}
						}
					}
					else {
						roff[y][x] = roff[y][x - 1];
					}
				}
				else if (y + max_abs_dy >= img->height) {
					*ptr++ = base;

					for (k = 0; k < inter_prd_order - 1; k++) {
						dy = idyx[k].y;
						dx = idyx[k].x;

						if(y + dy >= img->height || x + dx < 0 || x + dx >= img->width){
							*ptr++ = base;
						}
						else{
							*ptr++ = dy * img->width + dx + base;
						}
					}
				}
				else {
					roff[y][x] = roff[y - 1][x];
				}
			}
		}
	}

	return (roff);
}

/*---------------------------- init_encoder ------------------------*
 |  Function init_encoder
 |
 |  Purpose:  Initializes the encoder structure
 |
 |  Parameters:
 |      img		 		--> Video/Image file structure (IN)
 |		num_class		--> Number of classes to use (number of different predictors) (IN)
 |		num_group		--> Number of groups ???? (IN)
 |		prd_order		--> Order of the predictors (number of pixels to use) (IN)
 |		inter_prd_order	--> Order of the inter predictors (number of pixels to use) (IN)
 |		coef_precision	--> Precision of the coefficients (IN)
 |		f_huffman		--> Huffman coding flag (IN)
 |		quadtree_depth	--> Quadtree flag (IN)
 |		num_pmodel		--> Number of probability models (IN)
 |		pm_accuracy		--> Probability model accuracy (IN)
 |
 |  Returns:  ENCODER* --> returns a encoder type structure
 *-------------------------------------------------------------------*/
ENCODER *init_encoder(IMAGE *img, int num_class, int num_group, int prd_order, int inter_prd_order, int coef_precision, int f_huffman, int quadtree_depth, int num_pmodel, int pm_accuracy){
	//Declare new encoder struct
	ENCODER *enc;
	int x, y, i, j, f;
	double c;

	//Allocation of the memory for the encoder struct
	enc = (ENCODER *)alloc_mem(sizeof(ENCODER));

	//Copy of the video/image properties to encoder
	enc->height = img->height;
	enc->width = img->width;
	enc->frames = img->frames;
	enc->maxval = img->maxval;

	// Copy of the encoding parameters
	enc->num_class = malloc(enc->frames * sizeof(int));
	for(i = 0; i < enc->frames; i++){
		enc->num_class[i] = num_class; // M
	}
	enc->num_group = num_group; // ?? (= 16)
	enc->prd_order = prd_order; // K
	enc->inter_prd_order = inter_prd_order; // J
	enc->coef_precision = coef_precision; // P
	enc->f_huffman = f_huffman; // -h
	enc->quadtree_depth = quadtree_depth; // -f
	enc->num_pmodel = num_pmodel; // V
	enc->pm_accuracy = pm_accuracy; // A
	enc->maxprd = enc->maxval << enc->coef_precision; // ??

	// Alloc memory to predictors array
	enc->predictor = (int ***) alloc_mem(enc->frames * sizeof(int **));
	enc->predictor[0] = (int **)alloc_2d_array(enc->num_class[0], enc->prd_order, sizeof(int));
	for(f = 1; f < enc->frames; f++){
		enc->predictor[f] = (int **)alloc_2d_array(enc->num_class[f], enc->prd_order + enc->inter_prd_order, sizeof(int));
	}
	//enc->predictor = (int ***)alloc_3d_array(enc->num_class[0], enc->prd_order + enc->inter_prd_order, enc->frames, sizeof(int));

	// Alloc memory to ??
	enc->th = (int ***)alloc_3d_array(enc->num_class[0], enc->num_group, enc->frames, sizeof(int));

	// enc->th initializing cycle
	for(f = 0; f < enc->frames; f++){
		for (i = 0; i < enc->num_class[f]; i++){
			for (j = 0; j < enc->num_group - 1; j++){
				enc->th[f][i][j] = 0;
			}
			enc->th[f][i][enc->num_group - 1] = MAX_UPARA + 1;
		}
	}

	// More memory allocation
	enc->upara = (int ***)alloc_3d_array(enc->height, enc->width, enc->frames, sizeof(int));
	enc->prd = (int ***)alloc_3d_array(enc->height, enc->width, enc->frames, sizeof(int));

	if(enc->frames > 0){
		enc->roff = (int ****) alloc_mem(2 * sizeof(int ***));
		enc->roff[0] = init_ref_offset(img, enc->prd_order, 0);
		enc->roff[1] = init_ref_offset(img, enc->prd_order, enc->inter_prd_order);
	}
	else{
		enc->roff = (int ****) alloc_mem(1 * sizeof(int ***));
		enc->roff[0] = init_ref_offset(img, enc->prd_order, 0);
	}

	// Memory allocation for the original image
	enc->org = (int ***)alloc_3d_array(enc->height + 1, enc->width, enc->frames, sizeof(int));
	// Memory allocation for the error image ??
	enc->err = (int ***)alloc_3d_array(enc->height + 1, enc->width, enc->frames, sizeof(int));

	enc->qtctx = (int **) alloc_2d_array(enc->frames, QUADTREE_DEPTH << 3, sizeof(int));

	// Quadtree map what?
	if (enc->quadtree_depth > 0){
		enc->qtmap = (char ****)alloc_2d_array(enc->frames, QUADTREE_DEPTH, sizeof(char**));

		for(f = 0; f < enc->frames; f++){
			y = (enc->height + MAX_BSIZE - 1) / MAX_BSIZE;
			x = (enc->width + MAX_BSIZE - 1) / MAX_BSIZE;

			for (i = enc->quadtree_depth - 1; i >= 0; i--){
				enc->qtmap[f][i] = (char **)alloc_2d_array(y, x, sizeof(char));
				y <<= 1;
				x <<= 1;
			}
		}
	}

	// ?
	enc->ctx_weight = (int **) alloc_mem(enc->frames * sizeof(int *));
	for(f = 0; f < enc->frames; f++){
		enc->ctx_weight[f] = init_ctx_weight();
	}

	// Class and group arrays
	enc->class = (char ***)alloc_3d_array(enc->height, enc->width, enc->frames, sizeof(char));
	enc->group = (char ***)alloc_3d_array(enc->height, enc->width, enc->frames, sizeof(char));

	// Copy original images to encoder struct, initialize group array
	for(f = 0; f < enc->frames; f++){
		for (y = 0; y < enc->height; y++){
			for (x = 0; x < enc->width; x++){
				enc->group[f][y][x] = 0;
				enc->org[f][y][x] = img->val[f][y][x];
			}
		}

		// ?
		enc->org[f][enc->height][0] = (enc->maxval + 1) >> 1;
		enc->err[f][enc->height][0] = (enc->maxval + 1) >> 2;
	}

	// ?
	enc->uquant = (char ***)alloc_3d_array(enc->num_class[0], MAX_UPARA + 1, enc->frames, sizeof(char));

	for(f = 0; f < enc->frames; f++){
		for (i = 0; i < enc->num_class[0]; i++){
			for (j = 0; j <= MAX_UPARA; j++){
				enc->uquant[f][i][j] = enc->num_group - 1;
			}
		}
	}

	// ?
	enc->econv = (int ***)alloc_3d_array(enc->maxval+1, (enc->maxval<<1)+1, enc->frames, sizeof(int));
	enc->bconv = (img_t **)alloc_2d_array(enc->frames, enc->maxprd + 1, sizeof(img_t));
	enc->fconv = (img_t **)alloc_2d_array(enc->frames, enc->maxprd + 1, sizeof(img_t));

	enc->pmlist = (PMODEL ***) alloc_mem(enc->frames * sizeof(PMODEL **));
	for(f = 0; f < enc->frames; f++){
		enc->pmlist[f] = (PMODEL **) alloc_mem(enc->num_group * sizeof(PMODEL *));
	}

	enc->spm = (PMODEL *) alloc_mem(enc->frames * sizeof(PMODEL));

	for(f = 0; f < enc->frames; f++){
		enc->spm[f].freq = alloc_mem((MAX_SYMBOL * 2 + 1) * sizeof(uint));
		enc->spm[f].cumfreq = &(enc->spm[f].freq[MAX_SYMBOL]);
	}

	// Huffman coding
	if (enc->f_huffman == 1){
		enc->sigma = sigma_h;
	}
	else{
		enc->sigma = sigma_a;
	}

	// ?
	enc->mtfbuf = (int **) alloc_2d_array(enc->frames, enc->num_class[0], sizeof(int));
	//enc->coef_m = (int **) alloc_2d_array(enc->frames, enc->prd_order, sizeof(int));
	enc->coef_m = (int **) alloc_mem(enc->frames * sizeof(int *));
	enc->coef_m[0] = (int *) alloc_mem(enc->prd_order * sizeof(int));
	for (i = 0; i < enc->prd_order; i++){
		enc->coef_m[0][i] = 0;
	}
	for(f = 1; f < enc->frames; f++){
		enc->coef_m[f] = (int *) alloc_mem((enc->prd_order + enc->inter_prd_order) * sizeof(int));
		for (i = 0; i < enc->prd_order + enc->inter_prd_order; i++){
			enc->coef_m[f][i] = 0;
		}
	}

	enc->coef_cost = (cost_t ***) alloc_3d_array(16, MAX_COEF + 1, enc->frames, sizeof(cost_t));

	for(f = 0; f < enc->frames; f++){
		for (i = 0; i < 16; i++){
#ifdef OPT_SIDEINFO
			if (enc->f_huffman == 1){
				for (j = 0; j <= MAX_COEF; j++){
					enc->coef_cost[f][i][j] = ((j >> i) + i + 1);

					if (j > 0) enc->coef_cost[f][i][j] += 1.0;
				}
			}
			else{
				double p;
				set_spmodel(&enc->spm[f], MAX_COEF + 1, i);
				p = log(enc->spm[f].cumfreq[MAX_COEF + 1]);

				for (j = 0; j <= MAX_COEF; j++){
					enc->coef_cost[f][i][j] = (p - log(enc->spm[f].freq[j])) / log(2.0);

					if (j > 0) enc->coef_cost[f][i][j] += 1.0;
				}
			}
#else
			for (j = 0; j <= MAX_COEF; j++){
				enc->coef_cost[f][i][j] = 0;
			}
#endif
		}
	}

	enc->th_cost = (cost_t **) alloc_mem(enc->frames * sizeof(cost_t));
	for(f = 0; f < enc->frames; f++){
		enc->th_cost[f] = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
	}

	for(f = 0; f < enc->frames; f++){
		for (i = 0; i < MAX_UPARA + 2; i++){
			enc->th_cost[f][i] = 0;
		}
	}

	enc->class_cost = (cost_t **) alloc_mem(enc->frames * sizeof(cost_t));
	for(f = 0; f < enc->frames; f++){
		enc->class_cost[f] = (cost_t *)alloc_mem(enc->num_class[f] * sizeof(cost_t));
	}

	c = log((double)enc->num_class[0]) / log(2.0);

	for(f = 0; f < enc->frames; f++){
		for (i = 0; i < enc->num_class[f]; i++){
			enc->class_cost[f][i] = c;
		}
	}

	enc->qtflag_cost = (cost_t **) alloc_2d_array(enc->frames, QUADTREE_DEPTH << 3, sizeof(cost_t));

	for(f = 0; f < enc->frames; f++){
		for (i = 0; i < (QUADTREE_DEPTH << 3); i++){
			enc->qtflag_cost[f][i] = 1.0;
		}
	}

	enc->vlcs = (VLC ***) alloc_mem(enc->frames * sizeof(VLC **));
	enc->pmodels = (PMODEL ****) alloc_mem(enc->frames * sizeof(PMODEL ***));
	enc->rc = (RANGECODER **) alloc_mem(enc->frames * sizeof(RANGECODER *));
	enc->encval = (int ***) alloc_mem(enc->frames * sizeof(int **));

	return (enc);
}

/*------------------------------ init_class -------------------------*
 |  Function init_class
 |
 |  Purpose:  Sets the class number of each pixel in a given frame
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |		frame			--> Frame to be processed (IN)
 |
 |  Returns:  void
 *-------------------------------------------------------------------*/
void init_class(ENCODER *enc, int frame){
	int k, x, y, i, j, v, cl, sum, num_block;
	int *var, *tmp, **ptr;

	// Number of blocks in the frame
	num_block = enc->height * enc->width / (BASE_BSIZE * BASE_BSIZE);

	var = (int *)  alloc_mem(num_block * sizeof(int));
	ptr = (int **) alloc_mem(num_block * sizeof(int *));

	//Variance calculation for each block
	for (k = 0; k < num_block; k++){
		//Gives the position of the top right pixel of each block
		y = (k / (enc->width / BASE_BSIZE)) * BASE_BSIZE;
		x = (k % (enc->width / BASE_BSIZE)) * BASE_BSIZE;

		var[k] = 0;
		sum = 0;

		//Run each pixel in a block
		for (i = 0; i < BASE_BSIZE; i++){
			for (j = 0; j < BASE_BSIZE; j++){
				v = enc->org[frame][y + i][x + j];
				sum += v;
				var[k] += v * v;
			}
		}

		//Final result of the variance for one block
		var[k] -= sum * sum / (BASE_BSIZE * BASE_BSIZE);
		ptr[k] = &(var[k]);
	}

	//Variance sorting, from lowest to highest
	for (i = num_block - 1; i > 0; i--){
		for (j = 0; j < i; j++){
			if (*ptr[j] > *ptr[j + 1]){
				tmp = ptr[j];
				ptr[j] = ptr[j + 1];
				ptr[j + 1] = tmp;
			}
		}
	}


	for(k = 0; k < num_block; k++){
		//Calculates the number of the class for a block, given its position (1, 2, ...)
		cl = (k * enc->num_class[frame]) / num_block;
		//Determines the new position of a block after being sorted
		v = (int)(ptr[k] - var);

		//Gives the position of the top right pixel of each of the sorted blocks
		y = (v / (enc->width / BASE_BSIZE)) * BASE_BSIZE;
		x = (v % (enc->width / BASE_BSIZE)) * BASE_BSIZE;

		//Sets the class number for each pixel
		for(i = 0; i < BASE_BSIZE; i++){
			for(j = 0; j < BASE_BSIZE; j++){
				enc->class[frame][y + i][x + j] = cl;
			}
		}
	}

	free(ptr);
	free(var);
}

void set_cost_model(ENCODER *enc, int frame, int f_mmse){
	int gr, i, j, k;
	double a, b, c, var;
	PMODEL *pm;

	for(i = 0; i <= enc->maxval; i++){
		for(j = 0; j <= (enc->maxval << 1); j++){
			k = (i << 1) - j;
			enc->econv[frame][i][j] = (k > 0)? (k - 1) : (-k);
		}
	}

	enc->encval[frame] = enc->err[frame];

	for(gr = 0; gr < enc->num_group; gr++){
		var = enc->sigma[gr] * enc->sigma[gr];

		if(f_mmse){
			a = 0;
			b = 1.0;
		}
		else{
			a = 0.5 * log(2 * M_PI * var) / log(2.0);
			b = 1.0 / (2.0 * log(2.0) * var);
		}

		enc->pmlist[frame][gr] = pm = enc->pmodels[frame][gr][enc->num_pmodel >> 1];

		for(k = 0; k <= pm->size; k++){
			c = (double)k * 0.5 + 0.25;
			pm->cost[k] = a + b * c * c;
		}

		pm->subcost[0] = 0.0;
	}

	for (k = 0; k <= enc->maxprd; k++){
		enc->bconv[frame][k] = 0;
		enc->fconv[frame][k] = 0;
	}
	return;
}

void set_cost_rate(ENCODER *enc, int frame){
	int gr, k, i, j, mask, shift, num_spm;
	double a, c;
	PMODEL *pm;

	if(enc->pm_accuracy < 0){
		for(i = 0; i <= enc->maxval; i++){
			for (j = 0; j <= (enc->maxval << 1); j++){
				k = (j + 1) >> 1;
				enc->econv[frame][i][j] = e2E(i - k, k, j&1, enc->maxval);
			}
		}
	}

	if (enc->pm_accuracy < 0){
		num_spm = 1;
	}
	else{
		enc->encval[frame] = enc->org[frame];
		mask = (1 << enc->pm_accuracy) - 1;
		shift = enc->coef_precision - enc->pm_accuracy;

		for(k = 0; k <= enc->maxprd; k++){
			i = (enc->maxprd - k + (1 << shift) / 2) >> shift;
			enc->fconv[frame][k] = (i & mask);
			enc->bconv[frame][k] = (i >> enc->pm_accuracy);
		}

		num_spm = 1 << enc->pm_accuracy;
	}

	a = 1.0 / log(2.0);

	for(gr = 0; gr < enc->num_group; gr++){
		for(i = 0; i < enc->num_pmodel; i++){
			pm = enc->pmodels[frame][gr][i];

			if (enc->f_huffman == 1){
				for (k = 0; k < pm->size; k++) {
					pm->cost[k] = enc->vlcs[frame][gr][pm->id].len[k];
				}

				pm->subcost[0] = 0.0;
			}
			else if(enc->pm_accuracy < 0){
				for(k = 0; k < pm->size; k++){
					pm->cost[k] = -a * log(pm->freq[k]);
				}

				c = pm->cumfreq[enc->maxval + 1];
				pm->subcost[0] = a * log(c);
			}
			else{
				for(j = 0; j < num_spm; j++){
					for (k = 0; k < pm->size; k++){
						pm->cost[k] = -a * log(pm->freq[k]);
					}
					for (k = 0; k <= enc->maxval; k++){
						c = pm->cumfreq[k + enc->maxval + 1] - pm->cumfreq[k];
						pm->subcost[k] = a * log(c);
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
 |  Purpose: "Calculates the prediction of each pixel"
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |		frame			--> Frame to be processed (IN)
 |		tly				--> Top left Y boundary
 |		tlx				--> Top left X boundary
 |		bry				--> Bottom right Y boundary
 |		brx				--> Bottom right X boundary
 |
 |  Returns:  void
 *-------------------------------------------------------------------*/
void predict_region(ENCODER *enc, int frame, int tly, int tlx, int bry, int brx){
	int x, y, k, cl, prd, org;
	int *coef_p;
	int *prd_p;
	int *roff_p, **roff_pp;
	int *err_p, *org_p;
	char *class_p;

	int prd_order = enc->prd_order;
	int acesso = 0;

	if(frame > 0){
		prd_order += enc->inter_prd_order;
		acesso = 1;
	}

	//Runs all the pixels in a frame (due to the used boundaries)
	for (y = tly; y < bry; y++){
		class_p = &enc->class[frame][y][tlx];
		org_p = &enc->org[frame][y][tlx];
		roff_pp = &enc->roff[acesso][y][tlx];
		err_p = &enc->err[frame][y][tlx];
		prd_p = &enc->prd[frame][y][tlx];

		for (x = tlx; x < brx; x++){
			cl = *class_p++;
			roff_p = *roff_pp++;
			coef_p = enc->predictor[frame][cl];
			prd = 0;

			for (k = 0; k < prd_order; k++){
				prd += org_p[*roff_p++] * (*coef_p++);
			}

			org = *org_p++;
			*prd_p++ = prd;

			//Boundaries of the predicted values
			if (prd < 0) prd = 0;
			else if (prd > enc->maxprd) prd = enc->maxprd;

			prd >>= (enc->coef_precision - 1);
			*err_p++ = enc->econv[frame][org][prd];
		}
	}
}

int calc_uenc(ENCODER *enc, int frame, int y, int x){
	int u, k, *err_p, *roff_p, *wt_p;

	int acesso = 0;
	if(frame > 0){
		acesso = 1;
	}

	err_p = &enc->err[frame][y][x];
	roff_p = enc->roff[acesso][y][x];
	wt_p = enc->ctx_weight[frame];

	u = 0;

	for (k =0; k < NUM_UPELS; k++) {
		u += err_p[*roff_p++] * (*wt_p++);
	}

	u >>= 6;

	if (u > MAX_UPARA) u = MAX_UPARA;

	return (u);
}

cost_t calc_cost(ENCODER *enc, int frame, int tly, int tlx, int bry, int brx){
	cost_t cost;
	int x, y, u, cl, gr, prd, e, base, frac;
	int *upara_p, *prd_p, *encval_p;
	char *class_p, *group_p;
	PMODEL *pm;

	//    bry += (UPEL_DIST - 1);
	//    tlx -= (UPEL_DIST - 1);
	//    brx += (UPEL_DIST - 1);
	if (bry > enc->height) bry = enc->height;
	if (tlx < 0) tlx = 0;
	if (brx > enc->width) brx = enc->width;

	cost = 0;

	for(y = tly; y < bry; y++){
		class_p  = &enc->class[frame][y][tlx];
		group_p  = &enc->group[frame][y][tlx];
		upara_p  = &enc->upara[frame][y][tlx];
		encval_p = &enc->encval[frame][y][tlx];
		prd_p    = &enc->prd[frame][y][tlx];

		for (x = tlx; x < brx; x++){
			cl = *class_p++;
			*upara_p++ = u = calc_uenc(enc, frame, y, x);
			*group_p++ = gr = enc->uquant[frame][cl][u];
			e = *encval_p++;
			prd = *prd_p++;

			if (prd < 0) prd = 0;
			else if (prd > enc->maxprd) prd = enc->maxprd;

			base = enc->bconv[frame][prd];
			frac = enc->fconv[frame][prd];
			pm = enc->pmlist[frame][gr] + frac;
			cost += pm->cost[base + e] + pm->subcost[base];
		}
	}

	return (cost);
}

/*-------------------------- design_predictor -----------------------*
 |  Function design_predictor
 |
 |  Purpose:  Calculates the predictor coefficients for each class
 |
 |  Parameters:
 |		enc				--> Encoder structure (IN / OUT)
 |		frame			--> Frame to be processed (IN)
 |		f_mmse			--> Coding type (Huffman/Arithmetic) (IN)
 |
 |  Returns:  cost_t --> cost of the designed predictor
 *-------------------------------------------------------------------*/
cost_t design_predictor(ENCODER *enc, int frame, int f_mmse){
	double **mat, *weight, w, e, d, pivot;
	int x, y, i, j, k, cl, gr, pivpos, *index, *roff_p, *org_p;

	int prd_order = enc->prd_order;
	int acesso = 0;

	if(frame > 0){
		prd_order += enc->inter_prd_order;
		acesso = 1;
	}

	mat = (double **)alloc_2d_array(prd_order, prd_order + 1, sizeof(double));
	index = (int *)alloc_mem(sizeof(int) * prd_order);
	weight = (double *)alloc_mem(sizeof(double) * enc->num_group);

	//Weight choice
	for(gr = 0; gr < enc->num_group; gr++){
		if(f_mmse){
			weight[gr] = 1.0;
		}
		else{
			weight[gr] = 1.0 / (enc->sigma[gr] * enc->sigma[gr]);
		}
	}

	//Cycle that runs each class in a given frame
	for(cl = 0; cl < enc->num_class[frame]; cl++){
		//Variable mat initialization
		for(i = 0; i < prd_order; i++){
			for(j = 0; j <= prd_order; j++){
				mat[i][j] = 0.0;
			}
		}

		//Run each pixel
		for(y = 0; y < enc->height; y++){
			for(x = 0; x < enc->width; x++){
				if(enc->class[frame][y][x] != cl){ //Check if the pixel is member of the current class
					x += BASE_BSIZE - 1;
					continue;
				}

				gr = enc->group[frame][y][x];
				roff_p = enc->roff[acesso][y][x]; //This variable has the position of the reference pixels
				org_p = &enc->org[frame][y][x];

				//Fills the matrix mat with the reference pixels values
				for (i = 0; i < prd_order; i++){
					w = weight[gr] * org_p[roff_p[i]];

					for(j = i; j < prd_order; j++){
						mat[i][j] += w * org_p[roff_p[j]];
					}

					mat[i][prd_order] += w * org_p[0]; //Results column with the value of the pixel to be predicted
				}
			}
		}

		//Makes the bottom diagonal equal to the top diagonal
		for(i = 0; i < prd_order; i++){
			index[i] = i;

			for(j = 0; j < i; j++){
				mat[i][j] = mat[j][i];
			}
		}

		//"Gaussian elimination"
		for(i = 0; i < prd_order; i++){
			pivpos = i;
			pivot = fabs(mat[index[i]][i]);

			//Sorts the matrix pivots
			for(k = i + 1; k < prd_order; k++){
				if(fabs(mat[index[k]][i]) > pivot){
					pivot = fabs(mat[index[k]][i]);
					pivpos = k;
				}
			}

			//Change the indexes
			k = index[i];
			index[i] = index[pivpos];
			index[pivpos] = k;

			//If the pivot is not zero the actual Gaussian elimination is performed
			if(pivot > 1E-10){
				d = mat[index[i]][i];

				for(j = i; j <= prd_order; j++){
					mat[index[i]][j] /= d;
				}

				for(k = 0; k < prd_order; k++){
					if (k == i) continue;

					d = mat[index[k]][i];

					for(j = i; j <= prd_order; j++){
						mat[index[k]][j] -= d * mat[index[i]][j];
					}
				}
			}
		}

		w = (1 << enc->coef_precision);
		e = 0.0;

		//Rounds the coefficients and stores the error
		for(i = 0; i < prd_order; i++){
			if(fabs(mat[index[i]][i]) > 1E-10){ //Checks if a line is not zero
				d = mat[index[i]][prd_order] * w;
			}
			else{
				d = 0.0;
			}

			k = d;
			if (k > d) k--;

			if (k < -MAX_COEF){
				k = d = -MAX_COEF;
			}
			else if(k > MAX_COEF){
				k = d = MAX_COEF;
			}

			enc->predictor[frame][cl][i] = k; //Coefficient for a given reference in the predictor
			d -= k;
			e += d;
			mat[index[i]][prd_order] = d; //Stores the error
		}

		/* minimize mean rounding errors */
		k = e + 0.5;
		for (;k > 0; k--){
			d = 0;
			for(j = i = 0; i < prd_order; i++){
				if(mat[index[i]][prd_order] > d){
					d = mat[index[i]][prd_order];
					j = i;
				}
			}

			if (enc->predictor[frame][cl][j] < MAX_COEF) enc->predictor[frame][cl][j]++;

			mat[index[j]][prd_order] = 0;
		}
	}

	free(weight);
	free(index);
	free(mat);

	predict_region(enc, frame, 0, 0, enc->height, enc->width);

	return (calc_cost(enc, frame, 0, 0, enc->height, enc->width));
}

cost_t optimize_group(ENCODER *enc, int frame){
	cost_t cost, min_cost, **cbuf, *dpcost, *cbuf_p, *thc_p;
	int x, y, th1, th0, k, u, cl, gr, prd, e, base, frac;
	int **trellis, *tre_p;
	PMODEL *pm, **pm_p;

	trellis = (int **)alloc_2d_array(enc->num_group, MAX_UPARA + 2, sizeof(int));
	dpcost = (cost_t *)alloc_mem((MAX_UPARA + 2) * sizeof(cost_t));
	cbuf = (cost_t **)alloc_2d_array(enc->num_group, MAX_UPARA + 2, sizeof(cost_t));
	thc_p = enc->th_cost[frame];

	for (k = 0; k < MAX_UPARA + 2; k++) trellis[0][k] = 0;
	/* Dynamic programming */
	for (cl = 0; cl < enc->num_class[frame]; cl++) {
		for (gr = 0; gr < enc->num_group; gr++) {
			cbuf_p = cbuf[gr];

			for (u = 0; u < MAX_UPARA + 2; u++) {
				cbuf_p[u] = 0;
			}
		}

		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				if (enc->class[frame][y][x] == cl) {
					u = enc->upara[frame][y][x] + 1;
					e = enc->encval[frame][y][x];
					prd = enc->prd[frame][y][x];

					if (prd < 0) prd = 0;
					else if (prd > enc->maxprd) prd = enc->maxprd;

					base = enc->bconv[frame][prd];
					frac = enc->fconv[frame][prd];
					pm_p = enc->pmlist[frame];

					for (gr = 0; gr < enc->num_group; gr++) {
						pm = (*pm_p++) + frac;
						cbuf[gr][u] += pm->cost[base + e] + pm->subcost[base];
					}
				}
			}
		}

		for (gr = 0; gr < enc->num_group; gr++) {
			cbuf_p = cbuf[gr];

			for (u = 1; u < MAX_UPARA + 2; u++) {
				cbuf_p[u] += cbuf_p[u - 1];
			}
		}

		cbuf_p = cbuf[0];

		for (u = 0; u < MAX_UPARA + 2; u++) {
			dpcost[u] = cbuf_p[u];
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
			enc->th[frame][cl][gr - 1] = th1;
		}
	}

	/* set context quantizer */
	for (cl = 0; cl < enc->num_class[frame]; cl++) {
		u = 0;

		for (gr = 0; gr < enc->num_group; gr++) {
			for (; u < enc->th[frame][cl][gr]; u++) {
				enc->uquant[frame][cl][u] = gr;
			}
		}
	}

	/* renew groups */
	cost = 0;
	pm_p = enc->pmlist[frame];

	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			cl = enc->class[frame][y][x];
			u = enc->upara[frame][y][x];
			enc->group[frame][y][x] = gr = enc->uquant[frame][cl][u];
			e = enc->encval[frame][y][x];
			prd = enc->prd[frame][y][x];

			if (prd < 0) prd = 0;
			else if (prd > enc->maxprd) prd = enc->maxprd;

			base = enc->bconv[frame][prd];
			pm = pm_p[gr] + enc->fconv[frame][prd];
			cost += pm->cost[base + e] + pm->subcost[base];
		}
	}

	/* optimize probability models */
	if (enc->optimize_loop > 1 && enc->num_pmodel > 1) {
		if (enc->num_pmodel > MAX_UPARA + 2) {
			free(cbuf);
			cbuf = (cost_t **)alloc_2d_array(enc->num_group, enc->num_pmodel, sizeof(cost_t));
		}

		for (gr = 0; gr < enc->num_group; gr++) {
			for (k = 0; k < enc->num_pmodel; k++) {
				cbuf[gr][k] = 0;
			}
		}

		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				gr = enc->group[frame][y][x];
				e = enc->encval[frame][y][x];
				prd = enc->prd[frame][y][x];

				if (prd < 0) prd = 0;
				else if (prd > enc->maxprd) prd = enc->maxprd;

				base = enc->bconv[frame][prd];
				frac = enc->fconv[frame][prd];

				for (k = 0; k < enc->num_pmodel; k++) {
					pm = enc->pmodels[frame][gr][k] + frac;
					cbuf[gr][k] += pm->cost[base + e] + pm->subcost[base];
				}
			}
		}

		for (gr = 0; gr < enc->num_group; gr++) {
			pm = enc->pmodels[frame][gr][0];
			cost = cbuf[gr][0];

			for (k = 1; k < enc->num_pmodel; k++) {
				if (cost > cbuf[gr][k]) {
					cost = cbuf[gr][k];
					pm = enc->pmodels[frame][gr][k];
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

void set_prdbuf(ENCODER *enc, int frame, int **prdbuf, int **errbuf, int tly, int tlx, int bufsize){
	int x, y, brx, bry, cl, k, prd, *prdbuf_p, *errbuf_p, *coef_p;
	int buf_ptr, org, *org_p, *roff_p;

	int prd_order = enc->prd_order;
	int acesso = 0;

	if(frame > 0){
		prd_order += enc->inter_prd_order;
		acesso = 1;
	}

	brx = (tlx + bufsize < enc->width) ? (tlx + bufsize) : enc->width;
	bry = (tly + bufsize < enc->height) ? (tly + bufsize) : enc->height;

	for (cl = 0; cl < enc->num_class[frame]; cl++) {
		buf_ptr = bufsize * (tly % bufsize) + tlx % bufsize;

		for (y = tly; y < bry; y++) {
			prdbuf_p = &prdbuf[cl][buf_ptr];
			errbuf_p = &errbuf[cl][buf_ptr];
			buf_ptr += bufsize;
			org_p = &enc->org[frame][y][tlx];

			for(x = tlx; x < brx; x++){
				if(cl == enc->class[frame][y][x]){
					*prdbuf_p++ = enc->prd[frame][y][x];
					*errbuf_p++ = enc->err[frame][y][x];
					org_p++;
				}
				else {
					coef_p = enc->predictor[frame][cl];
					roff_p = enc->roff[acesso][y][x];
					prd = 0;

					for (k = 0; k < prd_order; k++) {
						prd += org_p[*roff_p++] * (*coef_p++);
					}

					org = *org_p++;
					*prdbuf_p++ = prd;

					if (prd < 0) prd = 0;
					else if (prd > enc->maxprd) prd = enc->maxprd;

					prd >>= (enc->coef_precision - 1);
					*errbuf_p++ = enc->econv[frame][org][prd];
				}
			}
		}
	}
}

int find_class(ENCODER *enc, int frame, int **prdbuf, int **errbuf, int tly, int tlx, int bry, int brx, int bufsize){
	cost_t cost, min_cost;
	int x, y, bufptr, cl, min_cl;
	char *class_p;
	int *prd_p, *prdbuf_p, *err_p, *errbuf_p;

	min_cost = 1E8;
	min_cl = 0;

	for (cl = 0; cl < enc->num_class[frame]; cl++) {
		bufptr = bufsize * (tly % bufsize) + tlx % bufsize;
		cost = enc->class_cost[frame][enc->mtfbuf[frame][cl]];

		for (y = tly; y < bry; y++) {
			class_p = &enc->class[frame][y][tlx];
			prd_p = &enc->prd[frame][y][tlx];
			prdbuf_p = &prdbuf[cl][bufptr];
			err_p = &enc->err[frame][y][tlx];
			errbuf_p = &errbuf[cl][bufptr];
			bufptr += bufsize;

			for (x = tlx; x < brx; x++) {
				*class_p++ = cl;
				*prd_p++ = *prdbuf_p++;
				*err_p++ = *errbuf_p++;
			}
		}

		cost += calc_cost(enc, frame, tly, tlx, bry, brx);

		if (cost < min_cost) {
			min_cost = cost;
			min_cl = cl;
		}
	}

	bufptr = bufsize * (tly % bufsize) + tlx % bufsize;
	for (y = tly; y < bry; y++) {
		class_p = &enc->class[frame][y][tlx];
		prd_p = &enc->prd[frame][y][tlx];
		prdbuf_p = &prdbuf[min_cl][bufptr];
		err_p = &enc->err[frame][y][tlx];
		errbuf_p = &errbuf[min_cl][bufptr];
		bufptr += bufsize;
		for (x = tlx; x < brx; x++) {
			*class_p++ = min_cl;
			*prd_p++ = *prdbuf_p++;
			*err_p++ = *errbuf_p++;
		}
	}

	return (min_cl);
}

cost_t vbs_class(ENCODER *enc, int frame, int **prdbuf, int **errbuf, int tly, int tlx, int blksize, int width, int level){
	int y, x, k, bry, brx, cl, bufsize, bufptr, ctx;
	int mtf_save[MAX_CLASS];
	char **qtmap;
	cost_t cost1, cost2, qtcost;
	char *class_p;
	int *prd_p, *prdbuf_p, *err_p, *errbuf_p;

	if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
		bufsize = MAX_BSIZE;
	}
	else {
		bufsize = BASE_BSIZE;
	}

	brx = (tlx + blksize < enc->width) ? (tlx + blksize) : enc->width;
	bry = (tly + blksize < enc->height) ? (tly + blksize) : enc->height;

	if (tlx >= brx || tly >= bry) return (0);

	for (k = 0; k < enc->num_class[frame]; k++) {
		mtf_save[k] = enc->mtfbuf[frame][k];
	}

	mtf_classlabel(enc->class[frame], enc->mtfbuf[frame], tly, tlx, blksize, width, enc->num_class[frame]);
	cl = find_class(enc, frame, prdbuf, errbuf, tly, tlx, bry, brx, bufsize);
	qtcost = enc->class_cost[frame][enc->mtfbuf[frame][cl]];
	if (level > 0) {
		/* context for quad-tree flag */
		ctx = 0;
		qtmap = enc->qtmap[frame][level - 1];
		y = (tly / MIN_BSIZE) >> level;
		x = (tlx / MIN_BSIZE) >> level;

		if (y > 0) {
			if (qtmap[y - 1][x] == 1) ctx++;
			if (brx < width && qtmap[y - 1][x + 1] == 1) ctx++;
		}

		if (x > 0 && qtmap[y][x - 1] == 1) ctx++;

		ctx = ((level - 1) * 4 + ctx) << 1;
		/* Quad-tree partitioning */
		cost1 = calc_cost(enc, frame, tly, tlx, bry, brx) + enc->class_cost[frame][enc->mtfbuf[frame][cl]] + enc->qtflag_cost[frame][ctx];
		blksize >>= 1;

		for (k = 0; k < enc->num_class[frame]; k++) {
			enc->mtfbuf[frame][k] = mtf_save[k];
		}

		qtcost = enc->qtflag_cost[frame][ctx + 1];
		qtcost += vbs_class(enc, frame, prdbuf, errbuf, tly, tlx, blksize, width, level - 1);
		qtcost += vbs_class(enc, frame, prdbuf, errbuf, tly, tlx+blksize, blksize, width, level - 1);
		qtcost += vbs_class(enc, frame, prdbuf, errbuf, tly+blksize, tlx, blksize, width, level - 1);
		qtcost += vbs_class(enc, frame, prdbuf, errbuf, tly+blksize, tlx+blksize, blksize, brx, level - 1);
		cost2 = calc_cost(enc, frame, tly, tlx, bry, brx) + qtcost;

		if (cost1 < cost2) {
			blksize <<= 1;

			for (k = 0; k < enc->num_class[frame]; k++) {
				enc->mtfbuf[frame][k] = mtf_save[k];
			}

			mtf_classlabel(enc->class[frame], enc->mtfbuf[frame], tly, tlx, blksize, width, enc->num_class[frame]);
			qtcost = enc->class_cost[frame][enc->mtfbuf[frame][cl]] + enc->qtflag_cost[frame][ctx];
			bufptr = bufsize * (tly % bufsize) + tlx % bufsize;

			for (y = tly; y < bry; y++) {
				class_p = &enc->class[frame][y][tlx];
				prd_p = &enc->prd[frame][y][tlx];
				prdbuf_p = &prdbuf[cl][bufptr];
				err_p = &enc->err[frame][y][tlx];
				errbuf_p = &errbuf[cl][bufptr];
				bufptr += bufsize;

				for (x = tlx; x < brx; x++) {
					*class_p++ = cl;
					*prd_p++ = *prdbuf_p++;
					*err_p++ = *errbuf_p++;
				}
			}

			tly = (tly / MIN_BSIZE) >> level;
			tlx = (tlx / MIN_BSIZE) >> level;
			bry = tly + 1;
			brx = tlx + 1;

			for (; level > 0; level--) {
				qtmap = enc->qtmap[frame][level - 1];
				for (y = tly; y < bry; y++) {
					for (x = tlx; x < brx; x++) {
						qtmap[y][x] = 0;
					}
				}

				tly <<= 1;
				tlx <<= 1;
				bry <<= 1;
				brx <<= 1;
			}
		}
		else {
			qtmap[y][x] = 1;
		}
	}
	return (qtcost);
}

cost_t optimize_class(ENCODER *enc, int frame){
	int y, x, i, blksize, level;
	int **prdbuf, **errbuf;

	if(enc->quadtree_depth >= 0 && enc->optimize_loop > 1){
		level = enc->quadtree_depth;
		blksize = MAX_BSIZE;
	}
	else{
		level = 0;
		blksize = BASE_BSIZE;
	}

	for(i = 0; i < enc->num_class[frame]; i++){
		enc->mtfbuf[frame][i] = i;
	}

	prdbuf =(int **) alloc_2d_array(enc->num_class[frame], blksize * blksize, sizeof(int));
	errbuf =(int **) alloc_2d_array(enc->num_class[frame], blksize * blksize, sizeof(int));

	for(y = 0; y < enc->height; y += blksize){
		for(x = 0; x < enc->width; x += blksize){
			set_prdbuf(enc, frame, prdbuf, errbuf, y, x, blksize);
			vbs_class(enc, frame, prdbuf, errbuf, y, x, blksize, enc->width, level);
		}
	}

	free(errbuf);
	free(prdbuf);

	return (calc_cost(enc, frame, 0, 0, enc->height, enc->width));
}

void optimize_coef(ENCODER *enc, int frame, int cl, int pos1, int pos2){
#define SEARCH_RANGE 11
#define SUBSEARCH_RANGE 3
	cost_t cbuf[SEARCH_RANGE * SUBSEARCH_RANGE], *cbuf_p;
	float *pmcost_p;
	int i, j, k, x, y, df1, df2, base;
	int prd, prd_f, shift, maxprd, *coef_p, *econv_p, *roff_p, *org_p;
	char *class_p;
	img_t *bconv_p, *fconv_p;
	PMODEL *pm, *pm_p;

	cbuf_p = cbuf;
	coef_p = enc->predictor[frame][cl];
	k = 0;

	int acesso = 0;
	if(frame > 0){
		acesso = 1;
	}

	for (i = 0; i < SEARCH_RANGE; i++) {
		y = coef_p[pos1] + i - (SEARCH_RANGE >> 1);

		if (y < 0) y = -y;

		if (y > MAX_COEF) y = MAX_COEF;

		for (j = 0; j < SUBSEARCH_RANGE; j++) {
			x = coef_p[pos2] - (i - (SEARCH_RANGE >> 1)) - (j - (SUBSEARCH_RANGE >> 1));

			if (x < 0) x = -x;

			if (x > MAX_COEF) x = MAX_COEF;

			cbuf_p[k++] = enc->coef_cost[frame][enc->coef_m[frame][pos1]][y] + enc->coef_cost[frame][enc->coef_m[frame][pos2]][x];
		}
	}

	bconv_p = enc->bconv[frame];
	fconv_p = enc->fconv[frame];
	maxprd = enc->maxprd;
	shift = enc->coef_precision - 1;

	for (y = 0; y < enc->height; y++) {
		class_p = enc->class[frame][y];

		for (x = 0; x < enc->width; x++) {
			if (cl != *class_p++) continue;

			roff_p = enc->roff[acesso][y][x];
			prd = enc->prd[frame][y][x];
			org_p = &enc->org[frame][y][x];
			pm_p = enc->pmlist[frame][(int)enc->group[frame][y][x]];
			df1 = org_p[roff_p[pos1]];
			df2 = org_p[roff_p[pos2]];
			prd_f = prd - (df1 - df2) * (SEARCH_RANGE >> 1) + df2 * (SUBSEARCH_RANGE >> 1);
			cbuf_p = cbuf;

			if (enc->pm_accuracy < 0) {
				econv_p = enc->econv[frame][*org_p];
				pmcost_p = pm_p->cost;

				for (i = 0; i < SEARCH_RANGE; i++) {
					for (j = 0; j < SUBSEARCH_RANGE; j++) {
						prd = prd_f;

						if (prd < 0) prd = 0;
						else if (prd > maxprd) prd = maxprd;

						(*cbuf_p++) += pmcost_p[econv_p[prd >> shift]];
						prd_f -= df2;
					}

					prd_f += df1 + df2 * (SUBSEARCH_RANGE - 1);
				}
			}
			else {
				for (i = 0; i < SEARCH_RANGE; i++) {
					for (j = 0; j < SUBSEARCH_RANGE; j++) {
						prd = prd_f;

						if (prd < 0) prd = 0;
						else if (prd > maxprd) prd = maxprd;

						base = bconv_p[prd];
						pm = pm_p + fconv_p[prd];
						(*cbuf_p++) += pm->cost[*org_p + base] + pm->subcost[base];
						prd_f -= df2;
					}

					prd_f += df1 + df2 * (SUBSEARCH_RANGE - 1);
				}
			}
		}
	}

	cbuf_p = cbuf;
	j = (SEARCH_RANGE * SUBSEARCH_RANGE) >> 1;

	for (i = 0; i < SEARCH_RANGE * SUBSEARCH_RANGE; i++) {
		if (cbuf_p[i] < cbuf_p[j]) {
			j = i;
		}
	}

	i = (j / SUBSEARCH_RANGE) - (SEARCH_RANGE >> 1);
	j = (j % SUBSEARCH_RANGE) - (SUBSEARCH_RANGE >> 1);
	y = coef_p[pos1] + i;
	x = coef_p[pos2] - i - j;

	if (y < -MAX_COEF) y = -MAX_COEF;
	else if (y > MAX_COEF) y = MAX_COEF;

	if (x < -MAX_COEF) x = -MAX_COEF;
	else if (x > MAX_COEF) x = MAX_COEF;

	i = y - coef_p[pos1];
	j = x - coef_p[pos2];

	if (i != 0 || j != 0) {
		for (y = 0; y < enc->height; y++) {
			class_p = enc->class[frame][y];

			for (x = 0; x < enc->width; x++) {
				if (cl == *class_p++) {
					roff_p = enc->roff[acesso][y][x];
					org_p = &enc->org[frame][y][x];
					enc->prd[frame][y][x] += org_p[roff_p[pos1]] * i + org_p[roff_p[pos2]] * j;
				}
			}
		}

		coef_p[pos1] += i;
		coef_p[pos2] += j;
	}
}

cost_t optimize_predictor(ENCODER *enc, int frame){
	int cl, k, pos1, pos2;
#ifndef RAND_MAX
#  define RAND_MAX 32767
#endif

	int prd_order = enc->prd_order;
	if(frame > 0) prd_order += enc->inter_prd_order;

	for (cl = 0; cl < enc->num_class[frame]; cl++) {
		for (k = 0; k < prd_order; k++) {
			retry:
			pos1 = (int)(((double)rand() * prd_order) / (RAND_MAX+1.0));
			pos2 = (int)(((double)rand() * prd_order) / (RAND_MAX+1.0));

			if (pos1 == pos2) goto retry;

			optimize_coef(enc, frame, cl, pos1, pos2);
		}
	}

	predict_region(enc, frame, 0, 0, enc->height, enc->width);
	return (calc_cost(enc, frame, 0, 0, enc->height, enc->width));
}

// Write bits to file
int putbits(FILE *fp, int n, uint x){
	static int bitpos = 8;
	static uint bitbuf = 0;
	int bits;

	bits = n;

	if (bits <= 0) return (0);

	while(n >= bitpos){
		n -= bitpos;

		if(n < 32){
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

void remove_emptyclass(ENCODER *enc, int frame){
	int cl, i, k, x, y;

	int prd_order = enc->prd_order;
	if(frame > 0) prd_order += enc->inter_prd_order;

	for (cl = 0; cl < enc->num_class[frame]; cl++) {
		enc->mtfbuf[frame][cl] = 0;
	}

	for (y = 0; y < enc->height; y += MIN_BSIZE) {
		for (x = 0; x < enc->width; x += MIN_BSIZE) {
			cl = enc->class[frame][y][x];
			enc->mtfbuf[frame][cl]++;
		}
	}

	for (i = cl = 0; i < enc->num_class[frame]; i++) {
		if (enc->mtfbuf[frame][i] == 0) {
			enc->mtfbuf[frame][i] = -1;
		}
		else {
			enc->mtfbuf[frame][i] = cl++;
		}
	}

	if (cl == enc->num_class[frame]) return;	/* no empty class */

	for (y = 0; y < enc->height; y++) {
		for (x = 0; x < enc->width; x++) {
			i = enc->class[frame][y][x];
			enc->class[frame][y][x] = enc->mtfbuf[frame][i];
		}
	}

	for (i = cl = 0; i < enc->num_class[frame]; i++) {
		if (enc->mtfbuf[frame][i] < 0) continue;
		if (cl != i) {
			for (k = 0; k < prd_order; k++) {
				enc->predictor[frame][cl][k] = enc->predictor[frame][i][k];
			}
			for (k = 0; k < enc->num_group - 1; k++) {
				enc->th[frame][cl][k] = enc->th[frame][i][k];
			}
		}
		cl++;
	}
	printf("M = %d\n", cl);
	enc->num_class[frame] = cl;
}

// Global header, only needs to be written once
int write_header(ENCODER *enc, FILE *fp){
	int bits, f;

	bits = putbits(fp, 16, MAGIC_NUMBER);
	bits += putbits(fp, 8, VERSION);
	bits += putbits(fp, 16, enc->width);
	bits += putbits(fp, 16, enc->height);
	bits += putbits(fp, 16, enc->maxval);
	bits += putbits(fp, 16, enc->frames);
	for(f = 0; f < enc->frames; f++){
		bits += putbits(fp, 8, enc->num_class[f]);
	}
	bits += putbits(fp, 4, 1);	/* number of components (1 = monochrome) */
	bits += putbits(fp, 6, enc->num_group);
	bits += putbits(fp, 7, enc->prd_order);
	bits += putbits(fp, 8, enc->inter_prd_order);
	bits += putbits(fp, 6, enc->num_pmodel - 1);
	bits += putbits(fp, 4, enc->coef_precision - 1);
	bits += putbits(fp, 3, enc->pm_accuracy + 1);
	bits += putbits(fp, 1, enc->f_huffman);
	bits += putbits(fp, 1, (enc->quadtree_depth < 0)? 0 : 1);

	return (bits);
}

int encode_golomb(FILE *fp, int m, int v){
	int bits, p;

	bits = p = (v >> m) + 1;

	while (p > 32) {
		putbits(fp, 32, 0);
		p -= 32;
	}

	putbits(fp, p, 1);	/* prefix code */
	putbits(fp, m, v);

	return (bits + m);
}

void set_qtindex(ENCODER *enc, int frame, int *index, uint *hist, int *numidx, int tly, int tlx, int blksize, int width, int level){
	int i, cl, x, y, ctx;
	char **qtmap;

	if (tly >= enc->height || tlx >= enc->width) return;

	if (level > 0) {
		/* context modeling for quad-tree flag */
		ctx = 0;
		qtmap = enc->qtmap[frame][level - 1];
		y = (tly / MIN_BSIZE) >> level;
		x = (tlx / MIN_BSIZE) >> level;

		if (y > 0) {
			if (qtmap[y - 1][x] == 1) ctx++;
			if (tlx + blksize < width && qtmap[y - 1][x + 1] == 1) ctx++;
		}

		if (x > 0 && qtmap[y][x - 1] == 1) ctx++;

		ctx = ((level - 1) * 4 + ctx) << 1;

		if (qtmap[y][x] == 1) {
			ctx++;
			index[(*numidx)++] = -(ctx + 1);
			enc->qtctx[frame][ctx]++;
			blksize >>= 1;
			set_qtindex(enc, frame, index, hist, numidx, tly, tlx, blksize, width, level - 1);
			set_qtindex(enc, frame, index, hist, numidx, tly, tlx + blksize, blksize, width, level - 1);
			set_qtindex(enc, frame, index, hist, numidx, tly + blksize, tlx, blksize, width, level - 1);
			width = tlx + blksize * 2;

			if (width >= enc->width) width = enc->width;

			set_qtindex(enc, frame, index, hist, numidx, tly + blksize, tlx + blksize, blksize, width, level - 1);

			return;
		}
		else {
			index[(*numidx)++] = -(ctx + 1);
			enc->qtctx[frame][ctx]++;
		}
	}

	cl = enc->class[frame][tly][tlx];
	mtf_classlabel(enc->class[frame], enc->mtfbuf[frame], tly, tlx, blksize, width, enc->num_class[frame]);
	i = enc->mtfbuf[frame][cl];
	index[(*numidx)++] = i;
	hist[i]++;

	return;
}

int encode_class(FILE *fp, ENCODER *enc, int frame){
	int i, j, k, numidx, blksize, level, x, y, ctx, bits, *index;
	uint *hist;
	cost_t cost;

#ifndef OPT_SIDEINFO
	if (fp == NULL) return(0);
#endif
	if (enc->quadtree_depth >= 0 && enc->optimize_loop > 1) {
		level = enc->quadtree_depth;
		blksize = MAX_BSIZE;
		numidx = 0;
		y = enc->height / MIN_BSIZE;
		x = enc->width / MIN_BSIZE;

		for (k = 0; k <= level; k++) {
			numidx += y * x;
			y = (y + 1) >> 1;
			x = (x + 1) >> 1;
		}

		for (k = 0; k < QUADTREE_DEPTH << 3; k++) {
			enc->qtctx[frame][k] = 0;
		}

	}
	else {
		level = 0;
		blksize = BASE_BSIZE;
		numidx = enc->height * enc->width / (BASE_BSIZE * BASE_BSIZE);
	}

	hist = (uint *)alloc_mem(enc->num_class[frame] * sizeof(uint));
	index = (int *)alloc_mem(numidx * sizeof(int));

	for (i = 0; i < enc->num_class[frame]; i++) {
		hist[i] = 0;
		enc->mtfbuf[frame][i] = i;
	}

	numidx = 0;

	for (y = 0; y < enc->height; y += blksize) {
		for (x = 0; x < enc->width; x += blksize) {
			set_qtindex(enc, frame, index, hist, &numidx, y, x, blksize, enc->width, level);
		}
	}

	bits = 0;

	if (enc->f_huffman == 1) {	/* Huffman */
		VLC *vlc;
		vlc = make_vlc(hist, enc->num_class[frame], 16);

		if (fp == NULL) {
			for (i = 0; i < enc->num_class[frame]; i++) {
				enc->class_cost[frame][i] = vlc->len[i];
			}

			for (j = 0; j < numidx; j++) {
				i = index[j];

				if (i < 0) {
					bits += 1;
				}
				else {
					bits += enc->class_cost[frame][i];
				}
			}
		}
		else {	/* actually encode */
			for (i = 0; i < enc->num_class[frame]; i++) {
				bits += putbits(fp, 4, vlc->len[i] - 1);
			}

			for (j = 0; j < numidx; j++) {
				i = index[j];

				if (i < 0) {
					ctx = -(i + 1);
					putbits(fp, 1, ctx & 1);
				}
				else {
					bits += putbits(fp, vlc->len[i], vlc->code[i]);
				}
			}
		}

		free_vlc(vlc);
	}
	else {			/* Arithmetic */
		PMODEL *pm;
		double p, c;
		int qtree_code[QUADTREE_DEPTH << 2], mtf_code[MAX_CLASS];

		/* context modeling for quad-tree flag */
		if (level > 0) {
			for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
				cost = INT_MAX;

				for (i = k = 0; i < 7; i++) {
					p = qtree_prob[i];
					c = -log(p) * (cost_t)enc->qtctx[frame][(ctx << 1) + 1] - log(1.0 - p) * (cost_t)enc->qtctx[frame][ctx << 1];

					if (c < cost) {
						k = i;
						cost = c;
					}
				}
				p = qtree_prob[qtree_code[ctx] = k];
				enc->qtflag_cost[frame][(ctx << 1) + 1] = -log(p) / log(2.0);
				enc->qtflag_cost[frame][ctx << 1] = -log(1.0 - p) / log(2.0);
			}
		}

		/* quantization of log-transformed probability */
		c = 0.0;

		for (i = 0; i < enc->num_class[frame]; i++) {
			c += (double)hist[i];
		}

		for (i = 0; i < enc->num_class[frame]; i++) {
			p = (double)hist[i] / c;

			if (p > 0.0) {
				mtf_code[i] = -log(p) / log(2.0) * (PMCLASS_LEVEL / PMCLASS_MAX);

				if (mtf_code[i] >= PMCLASS_LEVEL) {
					mtf_code[i] = PMCLASS_LEVEL - 1;
				}
			}
			else {
				mtf_code[i] = PMCLASS_LEVEL - 1;
			}

			p = exp(-log(2.0) * ((double)mtf_code[i] + 0.5) * PMCLASS_MAX / PMCLASS_LEVEL);
			enc->class_cost[frame][i] = -log(p) / log(2.0);
			hist[i] = p * (1 << 10);

			if (hist[i] <= 0) hist[i] = 1;
		}

		if (fp == NULL) {
			cost = 0.0;

			for (j = 0; j < numidx; j++) {
				i = index[j];

				if (i < 0) {
					ctx = -(i + 1);
					cost += enc->qtflag_cost[frame][ctx];
				}
				else {
					cost += enc->class_cost[frame][i];
				}
			}

			bits = (int)cost;
		}
		else {	/* actually encode */
			PMODEL cpm[1];
			/* additional info. */
			pm = &enc->spm[frame];

			if (level > 0) {
				set_spmodel(pm, 7, -1);

				for (ctx = 0; ctx < QUADTREE_DEPTH << 2; ctx++) {
					i = qtree_code[ctx];
					rc_encode(fp, enc->rc[frame], pm->cumfreq[i], pm->freq[i], pm->cumfreq[pm->size]);
				}
			}

			set_spmodel(pm, PMCLASS_LEVEL, -1);

			for (i = 0; i < enc->num_class[frame]; i++) {
				j = mtf_code[i];
				rc_encode(fp, enc->rc[frame], pm->cumfreq[j], pm->freq[j], pm->cumfreq[pm->size]);

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

			cpm->size = enc->num_class[frame];
			cpm->freq = (uint *)alloc_mem((cpm->size * 2 + 1) * sizeof(uint));
			cpm->cumfreq = &(cpm->freq[cpm->size]);
			cpm->cumfreq[0] = 0;

			for (i = 0; i < enc->num_class[frame]; i++) {
				cpm->freq[i] = hist[i];
				cpm->cumfreq[i + 1] = cpm->cumfreq[i] + cpm->freq[i];
			}

			for (j = 0; j < numidx; j++) {
				i = index[j];

				if (i < 0) {
					i = -(i + 1);
					ctx = i & (~1);
					rc_encode(fp, enc->rc[frame], pm->cumfreq[i] - pm->cumfreq[ctx], pm->freq[i], pm->cumfreq[ctx + 2] - pm->cumfreq[ctx]);
				}
				else {
					rc_encode(fp, enc->rc[frame], cpm->cumfreq[i], cpm->freq[i], cpm->cumfreq[cpm->size]);
				}
			}

			bits += enc->rc[frame]->code;
			enc->rc[frame]->code = 0;

			free(cpm->freq);
		}
	}

	free(index);
	free(hist);
	return (bits);
}

int encode_predictor(FILE *fp, ENCODER *enc, int frame){
	int cl, coef, sgn, k, m, min_m, bits;
	cost_t cost, min_cost, t_cost;

	int prd_order = enc->prd_order;
	if(frame > 0) prd_order += enc->inter_prd_order;

#ifndef OPT_SIDEINFO
	if (fp == NULL) return(0);
#endif
	t_cost = 0.0;

	for (k = 0; k < prd_order; k++) {
		min_cost = INT_MAX;

		for (m = min_m = 0; m < 16; m++) {
			cost = 0.0;

			for (cl = 0; cl < enc->num_class[frame]; cl++) {
				coef = enc->predictor[frame][cl][k];

				if (coef < 0) coef = -coef;

				cost += enc->coef_cost[frame][m][coef];
			}

			if (cost < min_cost) {
				min_cost = cost;
				min_m = m;
			}
		}

		t_cost += min_cost;
		enc->coef_m[frame][k] = min_m;
	}

	bits = t_cost;

	if (fp != NULL) {
		bits = 0;

		if (enc->f_huffman == 1) {	/* Huffman */
			for (k = 0; k < prd_order; k++) {
				bits += putbits(fp, 4, enc->coef_m[frame][k]);

				for (cl = 0; cl < enc->num_class[frame]; cl++) {
					coef = enc->predictor[frame][cl][k];
					sgn = (coef < 0)? 1 : 0;

					if (coef < 0) coef = -coef;

					bits += encode_golomb(fp, enc->coef_m[frame][k], coef);

					if (coef != 0) {
						bits += putbits(fp, 1, sgn);
					}
				}
			}
		} else {			/* Arithmetic */
			PMODEL *pm;
			pm = &enc->spm[frame];

			for (k = 0; k < prd_order; k++) {
				set_spmodel(pm, MAX_COEF + 1, enc->coef_m[frame][k]);
				rc_encode(fp, enc->rc[frame], enc->coef_m[frame][k], 1, 16);

				for (cl = 0; cl < enc->num_class[frame]; cl++) {
					coef = enc->predictor[frame][cl][k];
					sgn = (coef < 0)? 1 : 0;

					if (coef < 0) coef = -coef;

					rc_encode(fp, enc->rc[frame], pm->cumfreq[coef],  pm->freq[coef], pm->cumfreq[pm->size]);

					if (coef > 0) {
						rc_encode(fp, enc->rc[frame], sgn, 1, 2);
					}
				}
			}

			bits = enc->rc[frame]->code;
			enc->rc[frame]->code = 0;
		}
	}

	return (bits);
}

int encode_threshold(FILE *fp, ENCODER *enc, int frame){
	int cl, gr, i, k, m, min_m, bits;
	cost_t cost, min_cost;
	PMODEL *pm;

#ifndef OPT_SIDEINFO
	if (fp == NULL) return(0);
#endif
	if (enc->f_huffman == 1) {	/* Huffman */
		min_cost = INT_MAX;

		for (m = min_m = 0; m < 16; m++) {
			bits = 0;

			for (cl = 0; cl < enc->num_class[frame]; cl++) {
				k = 0;

				for (gr = 1; gr < enc->num_group; gr++) {
					i = enc->th[frame][cl][gr - 1] - k;
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
			enc->th_cost[frame][0] = 1.0;

			for (i = 1; i < MAX_UPARA + 2; i++) {
				enc->th_cost[frame][i] =((i - 1) >> min_m) + min_m + 1 + 1;
			}

			bits = min_cost;
		}
		else {
			bits = putbits(fp, 4, min_m);

			for (cl = 0; cl < enc->num_class[frame]; cl++) {
				k = 0;

				for (gr = 1; gr < enc->num_group; gr++) {
					i = enc->th[frame][cl][gr - 1] - k;

					if (i == 0) {
						bits += putbits(fp, 1, 0);
					} else {
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
					pm = enc->pmlist[frame][gr];
					bits += putbits(fp, k, pm->id);
				}
			}
		}
	}
	else {			/* Arithmetic */
		double p;
		pm = &enc->spm[frame];
		min_cost = INT_MAX;

		for (m = min_m = 0; m < 16; m++) {
			set_spmodel(pm, MAX_UPARA + 2, m);
			cost = 0.0;

			for (cl = 0; cl < enc->num_class[frame]; cl++) {
				k = 0;

				for (gr = 1; gr < enc->num_group; gr++) {
					i = enc->th[frame][cl][gr - 1] - k;
					p = (double)pm->freq[i] / (pm->cumfreq[pm->size - k]);
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
				enc->th_cost[frame][i] = (p - log(pm->freq[i])) / log(2.0);
			}

			bits = min_cost;
		}
		else {
			rc_encode(fp, enc->rc[frame], min_m, 1, 16);

			for (cl = 0; cl < enc->num_class[frame]; cl++) {
				k = 0;

				for (gr = 1; gr < enc->num_group; gr++) {
					i = enc->th[frame][cl][gr - 1] - k;
					rc_encode(fp, enc->rc[frame], pm->cumfreq[i],  pm->freq[i], pm->cumfreq[pm->size - k]);
					k += i;

					if (k > MAX_UPARA) break;
				}
			}

			if (enc->num_pmodel > 1) {
				for (gr = 0; gr < enc->num_group; gr++) {
					pm = enc->pmlist[frame][gr];
					rc_encode(fp, enc->rc[frame], pm->id, 1, enc->num_pmodel);
				}
			}

			bits = enc->rc[frame]->code;
			enc->rc[frame]->code = 0;
		}
	}

	return (bits);
}

int encode_image(FILE *fp, ENCODER *enc, int frame){
	int x, y, e, prd, base, bits, gr, cumbase;
	PMODEL *pm;

	bits = 0;
	if (enc->f_huffman == 1){	/* Huffman */
		VLC *vlc;

		for (y = 0; y < enc->height; y++){
			for (x = 0; x < enc->width; x++){
				gr = enc->group[frame][y][x];
				e = enc->encval[frame][y][x];
				pm = enc->pmlist[frame][gr];
				vlc = &enc->vlcs[frame][gr][pm->id];
				bits += putbits(fp, vlc->len[e], vlc->code[e]);
			}
		}
		
		if(frame == enc->frames - 1){
			putbits(fp, 7, 0);	/* flush remaining bits */
		}
	}
	else{			/* Arithmetic */
		for (y = 0; y < enc->height; y++){
			for (x = 0; x < enc->width; x++){
				gr = enc->group[frame][y][x];
				prd = enc->prd[frame][y][x];

				if (prd < 0) prd = 0;
				else if (prd > enc->maxprd) prd = enc->maxprd;

				e = enc->encval[frame][y][x];
				base = enc->bconv[frame][prd];
				pm = enc->pmlist[frame][gr] + enc->fconv[frame][prd];
				cumbase = pm->cumfreq[base];
				rc_encode(fp, enc->rc[frame], pm->cumfreq[base + e] - cumbase, pm->freq[base + e], pm->cumfreq[base + enc->maxval + 1] - cumbase);
			}
		}

		rc_finishenc(fp, enc->rc[frame]);
		bits += enc->rc[frame]->code;
	}

	return (bits);
}

void print_encoder(ENCODER *enc, int f){
	int i, j;
	printf("\n height = %d", enc->height);
	printf("\n width = %d", enc->width);
	printf("\n frames = %d", enc->frames);
	printf("\n maxval = %d", enc->maxval);
	printf("\n numclass = %d", enc->num_class[0]);
	printf("\n numgroup = %d", enc->num_group);
	printf("\n prd_prder = %d", enc->prd_order);
	printf("\n coefprecision = %d", enc->coef_precision);
	printf("\n num_pmodel = %d", enc->num_pmodel);
	printf("\n pmaccuracy = %d", enc->pm_accuracy);
	printf("\n maxprd = %d", enc->maxprd);
	printf("\n fhuffman = %d", enc->f_huffman);
	printf("\n quadtreedepth = %d", enc->quadtree_depth);
	printf("\n optimizeloop = %d", enc->optimize_loop);

	//	printf("\n predictor (%d x %d) =", enc->num_class, enc->prd_order);
	//	for(i = 0; i < enc->num_class - 1; i++)
	//		printf("\n");
	//		for(j = 0; j < enc->prd_order - 1; j++)
	//			printf("%d ", enc->predictor[f][i][j]);

	//	printf("\n th (%d x %d) =", enc->num_class, enc->prd_order);
	//		for(i = 0; i < enc->num_class - 1; i++)
	//			printf("\n");
	//			for(j = 0; j < enc->num_group - 1; j++)
	//				printf("%d ", enc->th[f][i][j]);

	//	printf("\n upara (%d x %d) =", enc->num_class, enc->prd_order);
	//	for(i = 0; i < enc->height - 1; i++)
	//		printf("\n");
	//		for(j = 0; j < enc->width - 1; j++)
	//			printf("%d ", enc->upara[f][i][j]);

	//	printf("\n prd (%d x %d) =", enc->num_class, enc->prd_order);
	//	for(i = 0; i < enc->height - 1; i++)
	//		printf("\n");
	//		for(j = 0; j < enc->width - 1; j++)
	//			printf("%d ", enc->prd[f][i][j]);

	//	printf("\n err (%d x %d) =", enc->num_class, enc->prd_order);
	//	for(i = 0; i < enc->height; i++)
	//		printf("\n");
	//		for(j = 0; j < enc->width - 1; j++)
	//			printf("%d ", enc->err[f][i][j]);

	//	printf("\n org =");
	//	for(i = 0; i < enc->height - 1; i++)
	//		printf("\n");
	//		for(j = 0; j < enc->width - 1; j++)
	//			printf("%d ", enc->org[f][i][j]);

	//	printf("\n ctxweight =");
	//	for(i = 0; i < NUM_UPELS - 1; i++)
	//		printf("%d ", enc->ctx_weight[f][i]);

	//	printf("\n qtctx =");
	//	for(i = 0; i < ((QUADTREE_DEPTH << 3) - 1); i++)
	//		printf("%d ", enc->qtctx[f][i]);

	printf("\n class =");
	for(i = 0; i < enc->height - 1; i++)
		printf("\n");
	for(j = 0; j < enc->width - 1; j++)
		printf("%d ", enc->class[f][i][j]);

	//	printf("\n group =");
	//	for(i = 0; i < enc->height - 1; i++)
	//		printf("\n");
	//		for(j = 0; j < enc->width - 1; j++)
	//			printf("%d ", enc->group[f][i][j]);

	//	printf("\n uquant =");
	//	for(i = 0; i < enc->num_class - 1; i++)
	//		printf("\n");
	//		for(j = 0; j < MAX_UPARA + 1 - 1; j++)
	//			printf("%d ", enc->uquant[f][i][j]);

	//	printf("\n econv =");
	//	for(i = 0; i < enc->maxval; i++)
	//		printf("\n");
	//		for(j = 0; j < (enc->maxval<<1); j++)
	//			printf("%d ", enc->econv[f][i][j]);

	//	printf("\n bconv =");
	//	for(i = 0; i < enc->maxprd; i++)
	//			printf("%d ", enc->bconv[f][i]);

	//	printf("\n fconv =");
	//	for(i = 0; i < enc->maxprd; i++)
	//			printf("%d ", enc->fconv[f][i]);

	//    PMODEL *spm;
	//    VLC ***vlcs;
	//    RANGECODER **rc;
	//    double *sigma;
	//    int **mtfbuf;
	//    int **coef_m;
	//    cost_t ***coef_cost;
	//    cost_t **th_cost;
	//    cost_t **class_cost;
	//    cost_t **qtflag_cost;//[QUADTREE_DEPTH << 3];
	printf("\n\n");
}

// remove_ext: removes the "extension" from a file spec.
//   mystr is the string to process.
//   dot is the extension separator.
//   sep is the path separator (0 means to ignore).
// Returns an allocated string identical to the original but
//   with the extension removed. It must be freed when you're
//   finished with it.
// If you pass in NULL or the new string can't be allocated,
//   it returns NULL.
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

// Main function
int main(int argc, char **argv){
	// Variable declaration
	cost_t cost, min_cost, side_cost;
	int i, j, f, k, x, y, cl, *bits, total_bits = 0, ***prd_save, ***th_save, resultado = 0;
	char ***class_save;
	double rate;
	IMAGE *img;
	ENCODER *enc;
	double elapse = 0.0;
	int f_mmse = 0;
	int f_optpred = 0;
	int f_huffman = 0;
	int quadtree_depth = QUADTREE_DEPTH;
	int num_class = NUM_CLASS;
	int num_group = NUM_GROUP;
	int prd_order = PRD_ORDER;
	int inter_prd_order = INTER_PRD_ORDER;
	int coef_precision = COEF_PRECISION;
	int num_pmodel = NUM_PMODEL;
	int pm_accuracy = PM_ACCURACY;
	int max_iteration = MAX_ITERATION;
	int frames = FRAMES;
	int height, width;
	char *infile, *outfile;
	char resfile[100];
	FILE *fp, *res;
	int predicao = 0;

	cpu_time();
	setbuf(stdout, 0);
	infile = outfile = NULL;

	// Read input parameters
	for(i = 1; i < argc; i++){
		if(argv[i][0] == '-'){
			switch(argv[i][1]){
			case 'H':
				height = atoi(argv[++i]);
				if(height <= 0){
					exit(-1);
				}
				break;
			case 'W':
				width = atoi(argv[++i]);
				if(width <= 0){
					exit(-2);
				}
				break;
			case 'F':
				frames = atoi(argv[++i]);
				if(frames <= 0){
					frames = FRAMES;
				}
				break;
			case 'p':
				predicao = atoi(argv[++i]);
				break;
			case 'M':
				num_class = atoi(argv[++i]);
				if(num_class <= 0 || num_class > 63){
					num_class = NUM_CLASS;
				}
				break;
			case 'K':
				prd_order = atoi(argv[++i]);
				if(prd_order <= 0 || prd_order > 72){
					prd_order = PRD_ORDER;
				}
				break;
			case 'J':
				inter_prd_order = atoi(argv[++i]);
//				if(prd_order <= 0 || prd_order > 72){
//					prd_order = PRD_ORDER;
//				}
				break;
			case 'P':
				coef_precision = atoi(argv[++i]);
				if(coef_precision <= 0 || coef_precision > 16){
					coef_precision = COEF_PRECISION;
				}
				break;
			case 'V':
				num_pmodel = atoi(argv[++i]);
				if(num_pmodel <= 0 || num_pmodel > 64){
					num_pmodel = NUM_PMODEL;
				}
				break;
			case 'A':
				pm_accuracy = atoi(argv[++i]);
				if(pm_accuracy < -1 || pm_accuracy > 6){
					pm_accuracy = PM_ACCURACY;
				}
				break;
			case 'I':
				max_iteration = atoi(argv[++i]);
				if(max_iteration <= 0){
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
			default:
				fprintf(stderr, "Unknown option: %s!\n", argv[i]);
				exit (1);
			}
		}
		else{
			if(infile == NULL){
				infile = argv[i];
			}
			else{
				outfile = argv[i];
			}
		}
	}

	//Memory allocation for the bit count
	bits = (int *) alloc_mem(frames * sizeof(int));

	//If the Huffman coding is used "turns off" the accuracy of the probabilities model
	if (f_huffman == 1) pm_accuracy = -1;
	//The accuracy of the probabilities model can't be greater than the coefficients precision
	if (pm_accuracy > coef_precision) pm_accuracy = coef_precision;

	// Help and version
	if(infile == NULL || outfile == NULL){
		printf(BANNER"\n", 0.1 * VERSION);
		printf("usage: encmrp [options] infile outfile\n");
		printf("options:\n");
		printf("    -H num  Height*\n");
		printf("    -W num  Width*\n");
		printf("    -F num  Frames*\n");
		printf("    -p num  Inter slice prediction type (ADICIONAR DESCRIÇÃO)\n");
		printf("    -M num  Number of predictors [%d]\n", num_class);
		printf("    -K num  Prediction order [%d]\n", prd_order);
		printf("    -J num  Inter Prediction order [%d]\n", inter_prd_order);
		printf("    -P num  Precision of prediction coefficients (fractional bits) [%d]\n", coef_precision);
		printf("    -V num  Number of probability models [%d]\n", num_pmodel);
		printf("    -A num  Accuracy of probability models [%d]\n", pm_accuracy);
		printf("    -I num  Maximum number of iterations [%d]\n", max_iteration);
		printf("    -m      Use MMSE predictors\n");
		printf("    -h      Use Huffman coding\n");
		printf("    -f      Fixed block-size for adaptive prediction\n");
		printf("    -o      Further optimization of predictors (experimental)\n");
		printf("infile:     Input file (must be in a raw PGM format)\n");
		printf("outfile:    Output file\n");
		printf("\nNote: * stands for a mandatory option.\n");
		exit(0);
	}

	// Read input file
	img = read_yuv(infile, height, width, frames);
	// Open output file
	fp = fileopen(outfile, "wb");

	// If the number of classes was not defined it is now
	k = img->width * img->height;
	if(num_class < 0){
		num_class = 10.4E-5 * k + 13.8;

		if (num_class > MAX_CLASS){
			num_class = MAX_CLASS;
		}
	}

	// If the prediction order was not defined it is now
	if(prd_order < 0){
		prd_order = 12.0E-5 * k + 17.2;

		for(i = 1; i < 8; i++){
			if(prd_order < (i + 1) * (i + 1)){
				prd_order = i * (i + 1);
				break;
			}
		}

		if (i >= 8) prd_order = 72;
	}

	printf("\nMRP-Video Encoder\n\n");
	// Print file characteristics to screen
	printf("%s -> %s (%dx%dx%d)\n", infile, outfile, img->width, img->height, img->frames);
	// Print coding parameters to screen
	printf("M = %d, K = %d, J = %d, P = %d, V = %d, A = %d\n", num_class, prd_order, inter_prd_order, coef_precision, num_pmodel, pm_accuracy);

	// Creates new ENCODER structure
	enc = init_encoder(img, num_class, num_group, prd_order, inter_prd_order, coef_precision, f_huffman, quadtree_depth, num_pmodel, pm_accuracy);
	
	//Prediction type selection
	if(predicao == 1){
		resultado = pixel_wise_predictor(enc);
	}
	else if(predicao == 2){
		resultado = gap_predictor(enc);
	}
	else if(predicao == 3){
		resultado = pixel_wise_search_predictor(enc, 3);
	}

	// Safeguard auxilliary variables
	prd_save = (int ***)alloc_mem(enc->frames * sizeof(int**));
	th_save = (int ***)alloc_mem(enc->frames * sizeof(int**));
	class_save = (char ***)alloc_mem(enc->frames * sizeof(char**));

	/* 1st loop */
//	#pragma omp parallel
//	{
//	#pragma omp for
	for(f = 0; f < enc->frames; f++){
		int prd_order = enc->prd_order;
		if(f > 0) prd_order += enc->inter_prd_order;

		// Initialize probability models
		enc->pmodels[f] = init_pmodels(enc->num_group, enc->num_pmodel, enc->pm_accuracy, NULL, enc->sigma, enc->maxval + 1);

		// Huffman coding
		if (enc->f_huffman == 1) {
			enc->vlcs[f] = init_vlcs(enc->pmodels[f], enc->num_group, enc->num_pmodel);
		}

		// Set cost model
		set_cost_model(enc, f, f_mmse);

		// Class initialization for each frame
		init_class(enc, f);

		prd_save[f] = (int **)alloc_2d_array(enc->num_class[f], prd_order, sizeof(int));
		th_save[f] = (int **)alloc_2d_array(enc->num_class[f], enc->num_group, sizeof(int));
		class_save[f] = (char **)alloc_2d_array(enc->height, enc->width, sizeof(char));

		//Number of the loop
		enc->optimize_loop = 1;
		min_cost = INT_MAX;

		for (i = j = 0; i < max_iteration; i++){
			cost = design_predictor(enc, f, f_mmse);
			cost = optimize_group(enc, f);
			cost = optimize_class(enc, f);

			if (cost < min_cost){
				min_cost = cost;
				j = i;

				for (y = 0; y < enc->height; y++){
					for (x = 0; x < enc->width; x++){
						class_save[f][y][x] = enc->class[f][y][x];
					}
				}

				for (cl = 0; cl < enc->num_class[f]; cl++){
					for (k= 0; k < prd_order; k++){
						prd_save[f][cl][k] = enc->predictor[f][cl][k];
					}

					for (k= 0; k < enc->num_group; k++){
						th_save[f][cl][k] = enc->th[f][cl][k];
					}
				}

			}

			if (i - j >= EXTRA_ITERATION) break;
			elapse += cpu_time();
		}

		for (y = 0; y < enc->height; y++) {
			for (x = 0; x < enc->width; x++) {
				enc->class[f][y][x] = class_save[f][y][x];
			}
		}

		for (cl = 0; cl < enc->num_class[f]; cl++) {
			for (k= 0; k < prd_order; k++) {
				enc->predictor[f][cl][k] = prd_save[f][cl][k];
			}

			for (k= 0; k < enc->num_group; k++) {
				enc->th[f][cl][k] = th_save[f][cl][k];
			}
		}

		set_cost_rate(enc, f);
		predict_region(enc, f, 0, 0, enc->height, enc->width);
		cost = calc_cost(enc, f, 0, 0, enc->height, enc->width);

		printf("Frame: %d --> Cost: %d\n", f, (int)cost);
	}
//	}

	/* 2nd loop */
	for(f = 0; f < enc->frames; f++){
		int prd_order = enc->prd_order;
		if(f > 0) prd_order += enc->inter_prd_order;

		enc->optimize_loop = 2;
		min_cost = INT_MAX;

		//printf("Frame: %d\n", f);
		for (i = j = 0; i < max_iteration; i++) {
			//printf("\t(%2d) cost =", i);

			if (f_optpred) {
				cost = optimize_predictor(enc, f);
				//printf(" %d ->", (int)cost);
			}

			side_cost = encode_predictor(NULL, enc, f);
			cost = optimize_group(enc, f);
			side_cost += encode_threshold(NULL, enc, f);
			//printf(" %d ->", (int)cost);
			cost = optimize_class(enc, f);
			side_cost += encode_class(NULL, enc, f);
			//printf(" %d (%d)", (int)cost, (int)side_cost);
			cost += side_cost;

			if (cost < min_cost) {
				//printf(" *\n");
				min_cost = cost;
				j = i;

				if (f_optpred) {
					for (y = 0; y < enc->height; y++) {
						for (x = 0; x < enc->width; x++) {
							class_save[f][y][x] = enc->class[f][y][x];
						}
					}
					for (cl = 0; cl < enc->num_class[f]; cl++) {
						for (k= 0; k < prd_order; k++) {
							prd_save[f][cl][k] = enc->predictor[f][cl][k];
						}
						for (k= 0; k < enc->num_group; k++) {
							th_save[f][cl][k] = enc->th[f][cl][k];
						}
					}
				}
			}
			else {
				//printf("\n");
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
			for (y = 0; y < enc->height; y++) {
				for (x = 0; x < enc->width; x++) {
					enc->class[f][y][x] = class_save[f][y][x];
				}
			}

			for (cl = 0; cl < enc->num_class[f]; cl++) {
				for (k= 0; k < prd_order; k++) {
					enc->predictor[f][cl][k] = prd_save[f][cl][k];
				}
				i = 0;

				for (k= 0; k < enc->num_group; k++) {
					enc->th[f][cl][k] = th_save[f][cl][k];
					for (; i < enc->th[f][cl][k]; i++) {
						enc->uquant[f][cl][i] = k;
					}
				}
			}
			predict_region(enc, f, 0, 0, enc->height, enc->width);
			calc_cost(enc, f, 0, 0, enc->height, enc->width);
			optimize_class(enc, f);
		}

		free(class_save[f]);
		free(prd_save[f]);
		free(th_save[f]);
		remove_emptyclass(enc, f);
		printf("Frame: %d --> Cost: %d (%d)\n", f, (int)cost, (int)side_cost);
	}

	free(class_save);
	free(prd_save);
	free(th_save);

	char *aux = remove_ext(outfile, '.', '/');
	sprintf(resfile, "res_%s.txt", aux);
	free(aux);

	res = fileopen(resfile, "w");
	fprintf(res, "MRP-Video version %d encoding results\n", VERSION);
	fprintf(res, "\tEncoded file: %s\n", infile);
	fprintf(res, "\tDimensions: %dx%dx%d\n", img->width, img->height, img->frames);
	fprintf(res, "---------------------------------------------\n");

	// Writes global header to file
	k = write_header(enc, fp);
	printf("------------------------------\n");
	printf("Header info.\t :%10d bits\n", k);
	fprintf(res, "Header info.\t :%10d bits\n", k);
	fprintf(res, "Frame\t\tBits\t\tBpp\n");
	total_bits = k;

	for(f = 0; f < enc->frames; f++){
		// Results output
		printf("------------------------------\n");
		printf("frame [%2d]\n", f);
		printf("header info.\t :%10d bits\n", k);

		if (enc->f_huffman == 0) {
			enc->rc[f] = rc_init();
			//putbits(fp, 7, 0);	/* byte alignment for the rangecoder */
		}

		bits[f] = 0;
		k = encode_class(fp, enc, f);
		bits[f] += k;
		printf("class info.\t :%10d bits\n", k);
		k = encode_predictor(fp, enc, f);
		bits[f] += k;
		printf("predictors\t :%10d bits\n", k);
		k = encode_threshold(fp, enc, f);
		bits[f] += k;
		printf("thresholds\t :%10d bits\n", k);
		k = encode_image(fp, enc, f);
		bits[f] += k;
		printf("pred. errors\t :%10d bits\n", k);
		printf("------------------------------\n");
		printf("total frame [%2d] :%10d bits\n", f, bits[f]);
		rate = (double)bits[f] / (enc->height * enc->width);
		printf("frame coding rate:%10.5f b/p\n", rate);
		fprintf(res, "%d\t%10d\t%10.5f\n", f, bits[f], rate);
	}

	for(f = 0; f < enc->frames; f++){
		total_bits += bits[f];
	}
	rate = (double)total_bits / (enc->height * enc->width * enc->frames);

	printf("------------------------------\n");
	printf("total\t\t :%10d bits\n", total_bits);
	printf("total coding rate:%10.5f b/p\n", rate);
	fclose(fp);
	elapse += cpu_time();
	printf("cpu time: %.2f sec.\n", elapse);

	fprintf(res, "---------------------------------------------\n");
	fprintf(res, "Total:\n\t%10d\t%10.5f\n", total_bits, rate);
	fprintf(res, "\nCPU time: %.2f sec.\n\n", elapse);

	if(resultado == 1){
		printf("The encoding process was lossy\n");
		fprintf(res, "The encoding process was lossy\n");
	}

	for(f = 0; f < enc->frames; f++){
		free(enc->ctx_weight[f]);

		for(i = enc->quadtree_depth - 1; i >= 0; i--){
			free(enc->qtmap[f][i]);
		}

		free(enc->pmlist[f]);
		free(enc->spm[f].freq);
		if(enc->f_huffman == 1){
			free(enc->vlcs[f]);
		}
		else{
			free(enc->rc[f]);
		}
		free(enc->th_cost[f]);
		free(enc->class_cost[f]);

		free(enc->pmodels[f]);
	}

	free(enc->roff[0]);
	if(enc->frames > 0) free(enc->roff[1]);

	free(img->val);
	free(img);
	free(bits);
	free(enc->num_class);
	free(enc->predictor);
	free(enc->th);
	free(enc->upara);
	free(enc->prd);
	free(enc->encval);
	free(enc->err);
	free(enc->org);
	free(enc->ctx_weight);
	free(enc->roff);
	free(enc->qtctx);
	free(enc->qtmap);
	free(enc->class);
	free(enc->group);
	free(enc->uquant);
	free(enc->econv);
	free(enc->bconv);
	free(enc->fconv);
	free(enc->pmodels);
	free(enc->pmlist);
	free(enc->spm);
	free(enc->vlcs);
	free(enc->rc);
	free(enc->mtfbuf);
	free(enc->coef_m);
	free(enc->coef_cost);
	free(enc->th_cost);
	free(enc->class_cost);
	free(enc->qtflag_cost);
	free(enc);
	fclose(res);
	return (0);
}
