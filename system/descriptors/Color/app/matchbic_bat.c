#include "colordescriptorslib.h"

int main(int argc, char** argv) {

    int n;			// Number of regions for each subimage
    Histogram **hist=NULL;	// List of histograms for each region
    char fvPath[100];
    FILE* output;

    ulong* distance;

    int r, i;

	if (argc != 3) {
		fprintf(stderr,"usage: match <curves_directory> <number_of_regions>\n");
		exit(-1);
	}

    n = atoi(argv[2]);
    hist = (Histogram**) calloc(n, sizeof(Histogram*));

    printf("Loading features...\n");
    for(r=0; r<n; r++) {
        sprintf(fvPath, "%s%d.fv", argv[1], r);
        hist[r] = ReadFileHistogram(fvPath);
        printf("%s read. \n", fvPath);
    }

    distance=(ulong*) calloc(n, sizeof(ulong));
    printf("Computing distances...\n");
    for(r=0; r<n; r++) {
        sprintf(fvPath, "%s%d.fv.distbin", argv[1], r);

        for(i=0; i<n; i++) {
            distance[i] = L1Distance(hist[r], hist[i]);
        }

        output=fopen(fvPath , "wb");
        fwrite(distance, sizeof(ulong), n, output);
        fclose(output);
    }

    for(r=0; r<n; r++) {
        DestroyHistogram(&hist[r]);
    }

    free(hist);
    free(distance);

    return(0);
}
