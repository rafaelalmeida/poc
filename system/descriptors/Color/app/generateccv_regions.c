#include "colordescriptorslib.h"

#define THRESHOLD 0.5

int isRelevant(Image *gt, Image *mask, int label) {
	
	int i,n;
	float rel, total;	
	
	n = gt->ncols * gt->nrows;
	
	rel = total = 0;
	i = 0;
	while(i<n) {
		if(mask->val[i] == label) {
			if(gt->val[i] > 0)
				rel++;
			total++;
		}
		i++;
	}	
	
	printf("%2.0f %2.0f\n", rel, total);
	
	if(rel/total >= THRESHOLD)
		return 1;
	return -1;
}

void WriteFeaturesByRegion(Histogram **hist, char* curve, char* groundtruth, Image *mask, int maxR) {
	
	Image *gt=NULL;
	int r, j;
	int hsize, label, value;
	
	FILE* output = fopen(curve, "w");
	
	gt = ReadImage(groundtruth);
	
	for(r = 0; r < maxR; r++) {
		
		label = isRelevant(gt,mask,r);		
		hsize = hist[r]->n;
		
		fprintf(output, "%d ", label);
		for(j = 0; j < hsize; j++) {
			value = hist[r]->v[j];
			if(value > 0)
				fprintf(output, "%d:%d ", j+1, value);
		}
		fprintf(output, "\n");
	}
	
	fclose(output);
	
	DestroyImage(&gt);	
}

void WriteFeaturesSVM(Histogram **hist, char* curve, char* groundtruth, Image *mask) {
	
	Image *gt=NULL;
	int i, j, n, aux;
	int hsize, label, value;
	
	FILE* output = fopen(curve, "w");
	
	gt = ReadImage(groundtruth);
	n = gt->ncols * gt->nrows;
	
	for(i = 0; i < n; i++) {
		
		label = 1;
		if(gt->val[i] == 0)
			label = -1;
		
		aux = mask->val[i];
		
		hsize = hist[mask->val[i]]->n;
		
		fprintf(output, "%d ", label);
		for(j = 0; j < hsize; j++) {
			value = hist[mask->val[i]]->v[j];
			if(value > 0)
				fprintf(output, "%d:%d ", j+1, value);
		}
		fprintf(output, "\n");
	}
	
	fclose(output);
	
	DestroyImage(&gt);
}

int main(int argc, char** argv)
{
  timer *tic, *toc;
  CImage *cimg=NULL;
  Image *mask=NULL;
  Image *rmask=NULL;
  Histogram **hist=NULL;
  int i, r, maxR;
  char outputfile[100];

  if (argc != 5) {
    fprintf(stderr,"usage: generateccv <image> <region mask> <curve> <ground truth>\n");
    exit(-1);
  }  

	cimg = ReadCImage(argv[1]);
	mask = ReadImage(argv[2]);

	maxR = MaximumValue(mask)+1;
	
	hist = (Histogram**) calloc(maxR, sizeof(Histogram*));
	
	for(r=0; r < maxR; r++) {
		
		hist[r] = NULL;
		
		rmask = CreateImage(mask->ncols,mask->nrows);

		for(i=0;i<(mask->ncols*mask->nrows);i++) {
			if(mask->val[i] == r)
				rmask->val[i] = 1;
		}
		  
		  //sprintf(outputfile, "%s%d.fv", argv[3], r);

		  tic = Tic();
		  hist[r] = CCV(cimg, rmask);
		  toc = Toc();
		  printf("\rCreateCCV for region %d in %f milliseconds\n",r, CTime(tic, toc));
		  //WriteFileHistogram(hist[r],outputfile);
		  
		  //sprintf(outputfile, "%d%s.pgm", r, argv[3]);
		  //WriteImage(rmask,outputfile);
		  
		  DestroyImage(&rmask);
	}

	printf("Writing features...\n");

	WriteFeaturesByRegion(hist, argv[3], argv[4], mask, maxR);

	for(r=0; r < maxR; r++)
		DestroyHistogram(&hist[r]);
	free(hist);

  //WriteFileHistogram(ccv,argv[2]);
  //DestroyHistogram(&ccv);
  DestroyCImage(&cimg);
  DestroyImage(&mask);

  return(0);
}
