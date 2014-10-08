/********************* Biblioteca de Projeto 2013 ***********************
* DESCRI��O: Biblioteca criada com fun��es comuns para a implementa��o  *
*            de v�rios programas.                                       *
*                                                                       *
* FUN��ES: le_yuv -> L� um ficheiro YUV para a mem�ria do programa      *
*          escreve_yuv -> Escreve um ficheiro YUV a partir das matrizes *
*                         fornecidas                                    *
*          mediana -> C�lcula a mediana do vetor fornecido              *
*          cmpfunc -> Retorna o resultado da compara��o de dois valores *
*          liberta_memoria -> liberta a mem�ria alocada para as         *
*                             componentes da imagem                     *
*                                                                       *
* AUTOR: Jo�o Santos                                                    *
* VERS�O: 3.0 - 28 de Abril de 2013                                     *
*************************************************************************/

#ifndef BIB_PROJ
#define BIB_PROJ

void le_yuv(int, int, int, char *, unsigned char ****, unsigned char ****, unsigned char ****);
void escreve_yuv(int, int, int, char *, unsigned char ***, unsigned char ***, unsigned char ***);
unsigned char mediana(int, unsigned char *);
int cmpfunc (const void*, const void*);
void liberta_memoria(int, int, int, unsigned char ***, unsigned char ***, unsigned char ***);

#endif
