#include "las.h"

void RGB2HSV_las(CImage *RGB, ImageHSV **HSV) {

    float r, g, b;
    float minVal, maxVal, delta;
    int i, j;

    for (i = 0; i < RGB->C[0]->nrows; i++) {
        for (j = 0; j < RGB->C[0]->ncols; j++) {

            //normaliza de 0 a 1
            r = (float) RGB->C[0]->val[j+RGB->C[0]->tbrow[i]]/255;
            g = (float) RGB->C[1]->val[j+RGB->C[1]->tbrow[i]]/255;
            b = (float) RGB->C[2]->val[j+RGB->C[2]->tbrow[i]]/255;

            minVal = min(r, g, b);
            maxVal = max(r, g, b);
            delta = maxVal - minVal;

            //H
            if (delta == 0) {
                HSV[i][j].H = 0;
            } else if (maxVal==r && g>=b) {
                HSV[i][j].H = 60*((g-b)/delta);
            } else if (maxVal==r && g<b) {
                HSV[i][j].H = 60*((g-b)/delta) + 360;
            } else if (maxVal==g) {
                HSV[i][j].H = 60*((b-r)/delta) + 120;
            } else if (maxVal==b) {
                HSV[i][j].H = 60* ((r-g)/delta) + 240;
            }

            //S
            if (maxVal==0) {
                HSV[i][j].S = 0;
            } else {
                HSV[i][j].S = delta/maxVal;
            }

            //V
            HSV[i][j].V = maxVal;

            //normalizando S e V entre 0 e 255
            HSV[i][j].S *= 255;
            HSV[i][j].V *= 255;

            //colocando H de 0 a 1
            //HSV[i][j].H /= 360;

        }
    }
}

void CalculaGradientes(Image *img, Image *g1, Image *g2, Image *g3, Image *g4) {

    int i,j;
    int va, vb; //usados para guardar os valores dos vizinhos (dependendo da direcao considerada)
                //facilita qdo sai dos limites da img
    int val_ij;

    //percorre imagem calculando g1, g2, g3 e g4
    for (i = 0; i < img->nrows; i++) {
        for (j = 0; j < img->ncols; j++) {

            val_ij = img->val[j+img->tbrow[i]]; //copia para uma variavel o valor do pixel atual
                                                //evita ter q ficar acessando a img toda hora

            //considera q pixels fora da imagem tem valor ZERO

            //g1  -- considera pixels da diagonal (descendo da esquerda pra direita)
            if (ValidPixel(img, j-1, i-1)) va=img->val[j-1+img->tbrow[i-1]];
            else va=0;

            if (ValidPixel(img, j+1, i+1)) vb=img->val[j+1+img->tbrow[i+1]];
            else vb=0;

            g1->val[j+g1->tbrow[i]] = abs(va - val_ij) + abs(val_ij - vb);


            //g2  -- considera pixels laterais
            if (ValidPixel(img, j-1, i)) va=img->val[j-1+img->tbrow[i]];
            else va=0;

            if (ValidPixel(img, j+1, i)) vb=img->val[j+1+img->tbrow[i]];
            else vb=0;

            g2->val[j+g2->tbrow[i]] = abs(va - val_ij) + abs(val_ij - vb);


            //g3  -- considera pixels da diagonal (subindo da esquerda pra direita)
            if (ValidPixel(img, j-1, i+1)) va=img->val[j-1+img->tbrow[i+1]];
            else va=0;

            if (ValidPixel(img, j+1, i-1)) vb=img->val[j+1+img->tbrow[i-1]];
            else vb=0;

            g3->val[j+g3->tbrow[i]] = abs(va - val_ij) + abs(val_ij - vb);


            //g4  -- considera pixels acima e abaixo
            if (ValidPixel(img, j, i+1)) va=img->val[j+img->tbrow[i+1]];
            else va=0;

            if (ValidPixel(img, j, i-1)) vb=img->val[j+img->tbrow[i-1]];
            else vb=0;

            g4->val[j+g4->tbrow[i]] = abs(va - val_ij) + abs(val_ij - vb);

            //if (g1->val[j+g1->tbrow[i]]>500 || g2->val[j+g2->tbrow[i]]>500 || g3->val[j+g3->tbrow[i]]>500 || g4->val[j+g4->tbrow[i]]>500)
            //    printf("g1=%d\tg2=%d\tg3=%d\tg4=%d\n", g1->val[j+g1->tbrow[i]], g2->val[j+g2->tbrow[i]], g3->val[j+g3->tbrow[i]], g4->val[j+g4->tbrow[i]]);
        }
    }
}

int Intervalo(int val) {
    //quantiza nao-uniformemente os 511 valores possiveis dos gradientes
    int lim1=25, lim2=50, lim3=75;

    if (val>=0 && val<lim1) {
        return 0;
    } else if (val>=lim1 && val<lim2) {
        return 1;
    } else if (val>=lim2 && val<lim3) {
        return 2;
    } else if (val>=lim3 && val<511) {
        return 3;
    } else {
        //printf("VALOR INVALIDO!");
        exit(1);
    }
    return -1;
}

float *HistogramaGradiente(Image *g1, Image *g2, Image *g3, Image *g4) {
    int i;
    float *hist=NULL;
    int bin=0; //bin a ser incrementado

    //aloca vetor histograma (calloc ja inicializa com zero)
    hist = (float*) calloc(QTD_BINS*QTD_BINS*QTD_BINS*QTD_BINS, sizeof(float));

    //percorre mapas de gradientes (todas tem o mesmo tamanho)
    for (i = 0; i < g1->nrows*g1->ncols; i++) {
        bin = Intervalo(g1->val[i]) +
              (Intervalo(g2->val[i])*QTD_BINS) +
              (Intervalo(g3->val[i])*QTD_BINS*QTD_BINS) +
              (Intervalo(g4->val[i])*QTD_BINS*QTD_BINS*QTD_BINS);

        hist[bin]++; //incrementa posicao correspondente
    }
    return hist;
}

void NormalizaHistograma(float *hist, int hist_size, int qtd_pixels) {
    int i;
    for (i=0; i<hist_size; i++) {
        hist[i] = hist[i]/qtd_pixels; //divide cada bin pelo tamanho da imagem
    }

}

/***************************************************************/

