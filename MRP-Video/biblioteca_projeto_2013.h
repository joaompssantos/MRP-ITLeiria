/********************* Biblioteca de Projeto 2013 ***********************
* DESCRIÇÃO: Biblioteca criada com funções comuns para a implementação  *
*            de vários programas.                                       *
*                                                                       *
* FUNÇÕES: le_yuv -> Lê um ficheiro YUV para a memória do programa      *
*          escreve_yuv -> Escreve um ficheiro YUV a partir das matrizes *
*                         fornecidas                                    *
*          mediana -> Cálcula a mediana do vetor fornecido              *
*          cmpfunc -> Retorna o resultado da comparação de dois valores *
*          liberta_memoria -> liberta a memória alocada para as         *
*                             componentes da imagem                     *
*                                                                       *
* AUTOR: João Santos                                                    *
* VERSÃO: 3.0 - 28 de Abril de 2013                                     *
*************************************************************************/

#ifndef BIB_PROJ
#define BIB_PROJ

void le_yuv(int, int, int, char *, unsigned char ****, unsigned char ****, unsigned char ****);
void escreve_yuv(int, int, int, char *, unsigned char ***, unsigned char ***, unsigned char ***);
unsigned char mediana(int, unsigned char *);
int cmpfunc (const void*, const void*);
void liberta_memoria(int, int, int, unsigned char ***, unsigned char ***, unsigned char ***);

#endif
