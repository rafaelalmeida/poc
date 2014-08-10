#include "qcch.h"

Image *RGB2Cinza(CImage *image)
{
    int i;
    Image *cinza=NULL;
    double aux=0;
    int img_size = (image->C[0]->ncols*image->C[0]->nrows);

    //cria imagem cinza
    cinza = CreateImage(image->C[0]->ncols,image->C[0]->nrows);

    //converte usando formula: 0.299*R + 0.587*G + 0.114*B
    for (i=0; i < img_size; i++) {
        aux = ( ((image->C[0]->val[i]) * 0.299)+((image->C[1]->val[i])*0.587) + ((image->C[2]->val[i])*0.114));
        cinza->val[i]= (int)aux;
    }

    return(cinza);
}

long int SomaJanela(int centroI, int centroJ, int raio, Image *img) {
    int k, m;
    long int y=0;

    for (k=-raio; k<=raio; k++) {
        for (m=-raio; m<=raio; m++) {
            if ( ValidPixel(img,(centroI+k),(centroJ+m)) ) { //verifica se o ponto eh valido
                y += img->val[(centroI+k)+img->tbrow[(centroJ+m)]];
            }
        }
    }
    //printf("SomaJanela[%d,%d]=%ld\t", centroI, centroJ, y);
    return y;
}

//image = cinza original
//t     = img resultado ja quantizada em 40 valores
//r     = raio considerado para o tamanho da janela
void Rate_Change_Gray(Image *image, Image *t, int r)
{
    int i,j,col=0,row=0;
    float hr, vr, dr, ar;  //directional rate change
    float vt=0.0;
    float y1, y2;

    col = (image->ncols);
    row = (image->nrows);

    for (j = r; j < row-r; j++) {
        for (i = r; i < col-r; i++) {

            //Calcula Hr - variacao horizontal
            y1 = (float) SomaJanela(i-r, j, r, image)/QTD_VIZ; //calcula media de cinza em  i-r,j
            y2 = (float) SomaJanela(i+r, j, r, image)/QTD_VIZ; //calcular media de cinza em i+r,j
            //usar media ao inves de valor do pixel
            //hr =  abs( (image->val[(i-r)+(image->tbrow[j])]) - (image->val[(i+r)+(image->tbrow[j])]) );
            hr =  fabsf(y1-y2);

            //Calcula Vr - variacao vertical
            y1 = (float) SomaJanela(i, j-r, r, image)/QTD_VIZ; //calcular media de cinza em i,j-r
            y2 = (float) SomaJanela(i, j+r, r, image)/QTD_VIZ; //calcular media de cinza em i,j+r
            //vr =  abs( (image->val[(i)+(image->tbrow[j-r])]) - (image->val[(i)+(image->tbrow[j+r])]) );
            vr = fabsf(y1-y2);

            //Calcula Dr - variacao diagonal
            y1 = (float) SomaJanela(i-r, j-r, r, image)/QTD_VIZ;  //calcular media de cinza em i-r,j-r
            y2 = (float) SomaJanela(i+r, j+r, r, image)/QTD_VIZ;  //calcular media de cinza em i+r,j+r
            //dr =  abs( (image->val[(i-r)+(image->tbrow[j-r])]) - (image->val[(i+r)+(image->tbrow[j+r])]) );
            dr = fabsf(y1-y2);

            //Calcula Ar - variacao anti-diagonal
            y1 = (float) SomaJanela(i+r, j-r, r, image)/QTD_VIZ;  //calcular media de cinza em i+r,j-r
            y2 = (float) SomaJanela(i-r, j+r, r, image)/QTD_VIZ;  //calcular media de cinza em i-r,j+r
            //ar =  abs( (image->val[(i+r)+(image->tbrow[j-r])]) - (image->val[(i-r)+(image->tbrow[j+r])]) );
            ar = fabsf(y1-y2);

            //vt eh a media dos valores
            vt = (hr + vr + dr + ar)/4.0;

            if (vt>255) printf("vt>255 -> %f\t hr=%f\t vr=%f\t dr=%f\t ar=%f\n", vt, hr, vr, dr, ar);

            //quantiza vt nao-uniformemente -- MUDAR ABAIXO!!!!!!!
               if( 0.0 < vt && vt <= 1.0 ){ (t->val[(i)+(t->tbrow[j])])=0; }
               else
               if( 1.0 < vt && vt <= 2.0 ){ (t->val[(i)+(t->tbrow[j])])=1; }
               else
               if( 2.0 < vt && vt <= 3.0 ){ (t->val[(i)+(t->tbrow[j])])=2; }
               else
               if( 3.0 < vt && vt <= 4.0 ){ (t->val[(i)+(t->tbrow[j])])=3; }
               else
               if( 4.0 < vt && vt <= 5.0 ){ (t->val[(i)+(t->tbrow[j])])=4; }
               else                       
               if( 5.0 < vt && vt <= 6.0 ){ (t->val[(i)+(t->tbrow[j])])=5; }
               else 
               if( 6.0 < vt && vt <= 7.0 ){ (t->val[(i)+(t->tbrow[j])])=6; }
               else
               if( 7.0 < vt && vt <= 8.0 ){ (t->val[(i)+(t->tbrow[j])])=7; }
               else
               if( 8.0 < vt && vt <= 9.0 ){ (t->val[(i)+(t->tbrow[j])])=8; }
               else
               if( 9.0 < vt && vt <= 10.0 ){ (t->val[(i)+(t->tbrow[j])])=9; }
               else            
               if( 10.0 < vt && vt <= 11.0 ){ (t->val[(i)+(t->tbrow[j])])=10; }
               else
               if( 11.0 < vt && vt <= 12.0 ){ (t->val[(i)+(t->tbrow[j])])=11; }
               else
               if( 12.0 < vt && vt <= 13.0 ){ (t->val[(i)+(t->tbrow[j])])=12; }
               else
               if( 13.0 < vt && vt <= 14.0 ){ (t->val[(i)+(t->tbrow[j])])=13; }
               else
               if( 14.0 < vt && vt <= 15.0 ){ (t->val[(i)+(t->tbrow[j])])=14; }
               else                       
               if( 15.0 < vt && vt <= 15.5 ){ (t->val[(i)+(t->tbrow[j])])=15; }
               else               
               if( 15.5 < vt && vt <= 17.5 ){ (t->val[(i)+(t->tbrow[j])])=16; }
               else
               if( 17.5 < vt && vt <= 19.5 ){ (t->val[(i)+(t->tbrow[j])])=17; }
               else
               if( 19.5 < vt && vt <= 21.5 ){ (t->val[(i)+(t->tbrow[j])])=18;  }
               else
               if( 21.5 < vt && vt <= 23.5 ){ (t->val[(i)+(t->tbrow[j])])=19; }
               else
               if( 23.5 < vt && vt <= 25.5 ){ (t->val[(i)+(t->tbrow[j])])=20; }
               else                       
               if( 25.5 < vt && vt <= 27.5 ){ (t->val[(i)+(t->tbrow[j])])=21; }
               else 
               if( 27.5 < vt && vt <= 29.5 ){ (t->val[(i)+(t->tbrow[j])])=22; }
               else
               if( 29.5 < vt && vt <= 31.5 ){ (t->val[(i)+(t->tbrow[j])])=23; }
               else
               if( 31.5 < vt && vt <= 33.5 ){ (t->val[(i)+(t->tbrow[j])])=24; }
               else
               if( 33.5 < vt && vt <= 35.5 ){ (t->val[(i)+(t->tbrow[j])])=25; }
               else               
               if( 35.5 < vt && vt <= 40.5 ){ (t->val[(i)+(t->tbrow[j])])=26; }
               else
               if( 40.5 < vt && vt <= 45.5 ){ (t->val[(i)+(t->tbrow[j])])=27; }
               else
               if( 45.5 < vt && vt <= 50.5 ){ (t->val[(i)+(t->tbrow[j])])=28; }
               else
               if( 50.5 < vt && vt <= 55.5 ){ (t->val[(i)+(t->tbrow[j])])=29; }
               else
               if( 55.5 < vt && vt <= 60.5 ){ (t->val[(i)+(t->tbrow[j])])=30; }
               else                       
               if( 60.5 < vt && vt <= 65.5 ){ (t->val[(i)+(t->tbrow[j])])=31; }
               else 
               if( 65.5 < vt && vt <= 70.5 ){ (t->val[(i)+(t->tbrow[j])])=32; }
               else
               if( 70.5 < vt && vt <= 75.5 ){ (t->val[(i)+(t->tbrow[j])])=33; }
               else
               if( 75.5 < vt && vt <= 80.5 ){ (t->val[(i)+(t->tbrow[j])])=34; }
               else
               if( 80.5 < vt && vt <= 85.5 ){ (t->val[(i)+(t->tbrow[j])])=35; }
               else 
               if( 85.5 < vt && vt <= 95.5 ){ (t->val[(i)+(t->tbrow[j])])=36; }
               else
               if( 95.5 < vt && vt <= 105.5 ){ (t->val[(i)+(t->tbrow[j])])=37; }
               else 
               if( 105.5 < vt && vt <= 115.5 ){ (t->val[(i)+(t->tbrow[j])])=38; }
               else
               if( 115.5 < vt && vt <= 255 ){ (t->val[(i)+(t->tbrow[j])])=39; }   
        }
    }  //fim - percorre img
}

//Calcula Normalized Quantized Compound Change Histogram
//retorn em h[] o histograma
void Change_Histogram(Image *t, double h[], int r)
{
    int i,tam,k,j;
    int aux;

    //ignora 2r linhas e 2r colunas
    tam = ((t->ncols)-2*r) * ((t->nrows)-2*r);

    //para cada um dos 40 valores,
    for(k=0; k<QTD_BINS; k++){ //h(t) onde t= 0...39
        aux = 0;
        //percorre a img t quantizada
        for (i = r; i < ((t->nrows)-r); i++) {
            for (j = r; j < ((t->ncols)-r); j++) {
                //calcula a funcao sigma:
                //  sigma(x)=1   se x=0,
                //  sigma(x)=0   caso contrario
                //no caso, x = t[i][j] - k
                if( ((t->val[(j)+(t->tbrow[i])]) - k ) == 0 ) {
                    aux+=1;
                }
            }
        }
        h[k] = (double)aux/tam;
    }
}

double *QCCH(CImage *cimg, int raio)
{
    Image *cinza=NULL;
    double *h=NULL; //vetor de caracteristicas

    Image *t = CreateImage(cimg->C[0]->ncols,cimg->C[0]->nrows);

    cinza = RGB2Cinza(cimg); //converte imagem de entrada

    h = (double*) calloc(QTD_BINS, sizeof(double));
    /*
    int i,j;
    printf("Imagem cinza\n");
    for (j = 0; j < t->nrows; j++) {
            for (i = 0; i < t->ncols; i++) {
                printf("%d\t", cinza->val[i+cinza->tbrow[j]]);
            }
            printf("\n");
    }*/

    Rate_Change_Gray(cinza, t, raio); //calcula imagem t com os valores de variacao de brilho nas vizinhancas
                                      //t ja retorna quantizada nao-uniformemente em 40 bins  (t->val[i] = [0,255])

    //Calcula Normalized Quantized Compound Change Histogram
    Change_Histogram(t, h, raio);

    //libera a memoria usada na funcao
    DestroyImage(&cinza);
    DestroyImage(&t);

    return h;
}

double Dist_Measure(double Ha[],double Hb[]) //array de QTD_BINS posicoes
{
    int i;
    double dab=0.0;

    for(i=0; i<QTD_BINS; i++) { //QTD_BINS bins
        dab+= fabs( Ha[i] - Hb[i] );
    }

    return dab;
}

/**********************************************************/
void ReadFile(char *filename, double h[])
{
    int i,n;
    FILE *fp;
    double c;

    fp = fopen(filename,"r");
    if (fp == NULL){
        fprintf(stderr,"Cannot open %s\n",filename);
        exit(-1);
    }
    fscanf(fp, "%d\n", &n);

    for (i=0; i < n; i++){
        fscanf(fp, "%lf\n", &c);
        h[i] = c; 
    }
    fclose(fp);
}

double *ReadFileBin(char *filename, int size) {
    FILE *fp;
    double *fv=NULL;

    fv = AllocDoubleArray(size);

    if ((fp = fopen(filename, "rb")) == NULL) {
        printf("ERRO CRIANDO ARQUIVO DE FV\n");
        exit(0);
    }
    fread(fv, sizeof(double), size, fp);
    fclose(fp);

    return fv;
}


void WriteFile(double h[], char *filename)
{
    FILE *fp;
    double *x;

    x=h;

    fp = fopen(filename,"w");
    if (fp == NULL){
        fprintf(stderr,"Cannot open %s\n",filename);
        exit(-1);
    }

    WriteStream(h,fp);
    fclose(fp);
}

void WriteStream(double *h, FILE *fp)
{
    int i;
    if (fp != NULL){
        fprintf(fp,"%d\n", QTD_BINS);
        for (i=0; i < QTD_BINS; i++) {
            fprintf(fp,"%f\n",h[i]);
        }
        fprintf(fp,"\n");
    }
}

void WriteFileBin(double *h, char *filename, int size) {
    FILE *fout;

    if ((fout = fopen(filename, "wb")) == NULL) {
        printf("ERRO CRIANDO ARQUIVO DE FV\n");
        exit(0);
    }
    fwrite(h, sizeof(double), size, fout);
    fclose(fout);
}
