#include "libcolordescriptors.h"

void Extraction(char *img_path, char *fv_path)
{
    double *h=NULL;

    CImage *cimg = ReadCImage(img_path); //leitura da imagem de entrada
    h = QCCH(cimg, RAIO); //calcula QCCH

    WriteFileBin(h, fv_path, QTD_BINS); //grava arquivo de FV

    //libera memoria
    DestroyCImage(&cimg);
    free(h);
}

void* LoadFV(char* fv_path) {
    return (void*) ReadFileBin(fv_path, QTD_BINS);
}

double Distance(void *fv1, void *fv2)
{
    double *ha = (double *) fv1, *hb = (double *) fv2;
    double distance;

    distance=Dist_Measure(ha,hb);
    //printf("%lf\n",distance);

    free(ha);
    free(hb);

    return distance;
}


