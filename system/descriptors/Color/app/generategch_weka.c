#include "colordescriptorslib.h"

#define THRESHOLD 0.8


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
	
	//printf("%2.0f %2.0f\n", rel, total);
	
	if(rel/total >= THRESHOLD)
		return 1;
	return 0;
}

int* computeRelevance(char* groundtruth, Image *mask) {
	
	int r, maxR;
	int* relevance;
	Image *gt=NULL;
	
	maxR = MaximumValue(mask)+1;
	relevance = (int*) calloc(maxR, sizeof(int));
	gt = ReadImage(groundtruth);
	
	for(r=0; r < maxR; r++)
	    relevance[r] = isRelevant(gt,mask,r);
	
	return relevance;
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
		
		for(j = 0; j < hsize; j++) {
			value = hist[r]->v[j];
			if(value > 0)
				fprintf(output, "%d ", value);
		}
		fprintf(output, "%d %d", r, label);
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
  Histogram *hist=NULL;
  int i, r, j, maxR, hsize, value, start;
  int *relevance;
  FILE* output;

  if (argc != 6) {
    fprintf(stderr,"usage: generategch <image> <region mask> <curve> <ground truth> <start_index>\n");
    exit(-1);
  }  

	cimg = ReadCImage(argv[1]);
	mask = ReadImage(argv[2]);
	start = atoi(argv[5]);
	
	maxR = MaximumValue(mask)+1;
	
	relevance = computeRelevance(argv[4], mask);
	rmask = CreateImage(mask->ncols,mask->nrows);
	
	for(r=start; r < maxR; r++) {
	  
		for(i=0;i<(mask->ncols*mask->nrows);i++) {
			if(mask->val[i] == r)
				rmask->val[i] = 1;
			else
				rmask->val[i] = 0;
		}
		  
		  //sprintf(outputfile, "%s%d.fv", argv[3], r);

		  tic = Tic();
		  hist = GCH(cimg, rmask);
		  toc = Toc();

		// writing
		hsize = hist->n;
		output = fopen(argv[3], "a");
		
		for(j = 0; j < hsize; j++) {
			value = hist->v[j];
//			if(value > 0)
				fprintf(output, "%d ", value);
		}
		fprintf(output, "%d %d\n", r, relevance[r]);
	
		printf("\rId: %4d   Concluido: %2.2f", r, (100.0 * r / (maxR-1)));
		
		fclose(output);
		DestroyHistogram(&hist);
	}
	DestroyImage(&rmask);
	
		
  //WriteFileHistogram(gch,argv[2]);
  //DestroyHistogram(&gch);
  DestroyCImage(&cimg);
  DestroyImage(&mask);

  return(0);
}
