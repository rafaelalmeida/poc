#include "libcolordescriptors.h"

int main(int argc, char** argv)
{
    double *h=NULL;

    if (argc != 3) {
        printf("\nNumero de entradas incorreto!");
        printf("\nUso correto: qcch_extraction <nome_arq_imagem_entrada> <nome_arq_imagem_saida> \n\n");
        return(1);
    }

    CImage *cimg = ReadCImage(argv[1]); //leitura da imagem de entrada
    h = QCCH(cimg, RAIO); //calcula QCCH

    WriteFileBin(h, argv[2], QTD_BINS); //grava arquivo de FV

    //libera memoria
    DestroyCImage(&cimg);
    free(h);

    return 0;
}


