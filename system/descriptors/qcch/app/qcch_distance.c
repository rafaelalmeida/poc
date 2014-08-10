#include "libcolordescriptors.h"

int main(int argc, char** argv)
{
    double *ha=NULL, *hb=NULL;
    double distance;

    if (argc != 3) {
        fprintf(stderr,"usage: qcch_distance <fv1_path> <fv2_path>\n");
        exit(-1);
    }

    ha = ReadFileBin(argv[1], QTD_BINS);
    hb = ReadFileBin(argv[2], QTD_BINS);

/*
    int i;
    printf("hist1=\n");
    for (i=0; i<QTD_BINS; i++) {
        printf("h1[%d]=%lf\t", i, ha[i]);
    }
    printf("\n");

    printf("hist2=\n");
    for (i=0; i<QTD_BINS; i++) {
        printf("h2[%d]=%lf\t", i, hb[i]);
    }
    printf("\n");
*/
    distance=Dist_Measure(ha,hb);

    printf("%lf\n",distance);

    free(ha);
    free(hb);

    return(0);
}
