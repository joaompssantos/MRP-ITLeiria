/********************* Biblioteca de Projeto 2013 ***********************
* DESCRIÇÃO: Biblioteca criada com funções comuns para a implementação  *
*            de vários programas.                                       *
*                                                                       *
* FUNÇÕES: le_yuv -> Lê um ficheiro YUV para a memória do programa      *
*          escreve_yuv -> Escreve um ficheiro YUV a partir das matrizes *
*                         fornecidas                                    *
*                                                                       *
* AUTOR: João Santos                                                    *
* VERSÃO: 3.0 - 28 de Abril de 2013                                     *
*************************************************************************/

#include "biblioteca_projeto_2013.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/***************************** le_yuv() *********************************
* DESCRIÇÃO: Função que lê uma imagem YUV para matrizes do programa     *
*                                                                       *
* INPUT: altura --> Altura em pixéis da imagem a ler                    *
*        largura --> Largura em pixéis da imagem a ler                  *
*        frames --> Número de frames do ficheiro a ler                  *
*        *ler --> Nome ou caminho da imagem a ler                       *
*        ****matrizY --> Matriz dos valores da componente Y da imagem   *
*        ****matrizU --> Matriz dos valores da componente U da imagem   *
*        ****matrizV --> Matriz dos valores da componente V da imagem   *
* RETURN: void                                                          *
*                                                                       *
* AUTOR: João Santos                                                    *
* VERSÃO: 3.0 - 12 de Maio de 2013                                      *
*************************************************************************/
void le_yuv(int altura, int largura, int frames, char* ler, unsigned char**** matrizY, unsigned char**** matrizU, unsigned char**** matrizV){
    FILE *le;

    int conta = 0, n = 0;

    unsigned char ***mY = (unsigned char ***) calloc(frames, sizeof(unsigned char **));
    for(n = 0; n < frames; n++){
        mY[n] = (unsigned char **) calloc(altura, sizeof(unsigned char *));
        for(conta = 0; conta < altura; conta++){
            mY[n][conta] = (unsigned char *) calloc(largura, sizeof(unsigned char));
        }
    }

    unsigned char ***mU = (unsigned char ***) calloc(frames, sizeof(unsigned char **));
    for(n = 0; n < frames; n++){
        mU[n] = (unsigned char **) calloc(altura / 2, sizeof(unsigned char *));
        for(conta = 0; conta < altura / 2; conta++){
            mU[n][conta] = (unsigned char *) calloc(largura / 2, sizeof(unsigned char));
        }
    }

    unsigned char ***mV = (unsigned char ***) calloc(frames, sizeof(unsigned char **));
    for(n = 0; n < frames; n++){
        mV[n] = (unsigned char **) calloc(altura / 2, sizeof(unsigned char *));
        for(conta = 0; conta < altura / 2; conta++){
            mV[n][conta] = (unsigned char *) calloc(largura / 2, sizeof(unsigned char));
        }
    }

    le = fopen(ler, "rb");

    if (le != NULL){
        for(n = 0; n < frames; n++){
            for(conta = 0; conta < altura; conta++){
                if(fread(mY[n][conta], sizeof(unsigned char), largura, le) != largura){
                    printf("\n Ocorreu um erro na leitura do ficheiro.");
                    exit(1);
                }
            }

            for(conta = 0; conta < altura / 2; conta++){
                if(fread(mU[n][conta], sizeof(unsigned char), largura / 2, le) != largura / 2){
                    printf("\n Ocorreu um erro na leitura do ficheiro.");
                    exit(2);
                }
            }

            for(conta = 0; conta < altura / 2; conta++){
                if(fread(mV[n][conta], sizeof(unsigned char), largura / 2, le) != largura / 2){
                    printf("\n Ocorreu um erro na leitura do ficheiro.");
                    exit(3);
                }
            }
        }
        fclose(le);
    }
    else{
        printf("\n Nao foi possivel ler o ficheiro. O programa vai ser terminado.");
        exit(-2);
    }

    *matrizY = mY;
    *matrizU = mU;
    *matrizV = mV;
}

/*************************** escreve_yuv() ******************************
* DESCRIÇÃO: Função que escreve um ficheiro YUV a partir das matrizes   *
*            do programa                                                *
*                                                                       *
* INPUT: altura --> Altura em pixéis da imagem a escrever               *
*        largura --> Largura em pixéis da imagem a escrever             *
*        frames --> Número de frames do ficheiro a escrever             *
*        *nome --> Nome ou caminho da imagem a escrever                 *
*        ***matrizY --> Matriz dos valores da componente Y da imagem    *
*        ***matrizU --> Matriz dos valores da componente U da imagem    *
*        ***matrizV --> Matriz dos valores da componente V da imagem    *
* RETURN: void                                                          *
*                                                                       *
* AUTOR: João Santos                                                    *
* VERSÃO: 3.0 - 12 de Maio de 2013                                      *
*************************************************************************/
void escreve_yuv(int altura, int largura, int frames, char *nome, unsigned char ***matrizY, unsigned char ***matrizU, unsigned char ***matrizV){
    FILE *escreve;
    int i = 0, n = 0;

    escreve = fopen(nome, "wb");

    if (escreve != NULL){
        for(n = 0; n < frames; n++){
            for(i = 0; i < altura; i++){
                fwrite(matrizY[n][i], sizeof(unsigned char), largura, escreve);
            }
            for(i = 0; i < altura / 2; i++){
                fwrite(matrizU[n][i], sizeof(unsigned char), largura / 2, escreve);
            }
            for(i = 0; i < altura / 2; i++){
                fwrite(matrizV[n][i], sizeof(unsigned char), largura / 2, escreve);
            }
        }
    }
    else{
        printf("\n Nao foi possivel escrever a imagem. O programa vai ser terminado.");
        exit(-3);
    }

    fclose(escreve);
}

/***************************** mediana() ********************************
* DESCRIÇÃO: Função que escreve um ficheiro YUV a partir das matrizes   *
*            do programa                                                *
*                                                                       *
* INPUT: tamanho --> Número de elementos do vetor do qual se calcula a  *
*                    mediana                                            *
*        *janela --> Vetor do qual se calcula a mediana                 *
*                                                                       *
* RETURN: median --> Resultado da mediana calculada                     *
*                                                                       *
* AUTOR: João Santos                                                    *
* VERSÃO: 1.0 - 21 de Abril de 2013                                     *
*************************************************************************/
unsigned char mediana(int tamanho, unsigned char *janela){
    unsigned char median = 0;

    //Ordenação crescente da janela
    qsort(janela, tamanho, sizeof(unsigned char), cmpfunc);
    //Cálculo da mediana
    if(tamanho % 2 == 0){
        median = ((janela[(tamanho / 2) - 1] + janela[tamanho / 2]) / 2) + 0.5f;
    }
    else{
        median = janela[tamanho / 2];
    }

    return median;
}

/***************************** cmpfunc() ********************************
* DESCRIÇÃO: Função que compara dois valores                            *
*                                                                       *
* AUTOR:                                                                *
* http://www.tutorialspoint.com/c_standard_library/c_function_qsort.htm *
*                                                                       *
* ACEDIDO A: 18 de Abril de 2013                                        *
*************************************************************************/
int cmpfunc (const void * a, const void * b){
   return ( *(unsigned char*)a - *(unsigned char*)b );
}

/***************************** mediana() ********************************
* DESCRIÇÃO: Função que desaloca as matrizes fornecidas                 *
*                                                                       *
* INPUT: altura --> Altura em pixéis da matriz Y                        *
*        largura --> Largura em pixéis da imagem matriz Y               *
*        frames --> Número de frames das matrizes                       *
*        ****matrizY --> Endereço da matriz da componente Y da imagem   *
*        ****matrizU --> Endereço da matriz da componente U da imagem   *
*        ****matrizV --> Endereço da matriz da componente V da imagem   *
*                                                                       *
* RETURN: void                                                          *
*                                                                       *
* AUTOR: João Santos                                                    *
* VERSÃO: 1.0 - 28 de Abril de 2013                                     *
*************************************************************************/
void liberta_memoria(int altura, int largura, int frames, unsigned char ***matrizY, unsigned char ***matrizU, unsigned char ***matrizV){
    int i = 0, j = 0;

    //Ciclos de desalocação de memória
    if(matrizY != NULL){
        for(i = 0; i < frames; i++){
            for(j = 0; j < altura; j++){
                free(matrizY[i][j]);
            }
            free(matrizY[i]);
        }
        free(matrizY);
    }

    if(matrizU != NULL){
        for(i = 0; i < frames; i++){
            for(j = 0; j < altura / 2; j++){
                free(matrizU[i][j]);
            }
            free(matrizU[i]);
        }
        free(matrizU);
    }

    if(matrizV != NULL){
        for(i = 0; i < frames; i++){
            for(j = 0; j < altura / 2; j++){
                free(matrizV[i][j]);
            }
            free(matrizV[i]);
        }
        free(matrizV);
    }
}
