#include "colordescriptorslib.h"

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
  int i, r, maxR;
  char outputfile[100];

  if (argc != 4) {
    fprintf(stderr,"usage: generatebic <image> <region mask> <curve>\n");
    exit(-1);
  }  

	cimg = ReadCImage(argv[1]);
	mask = ReadImage(argv[2]);

	maxR = MaximumValue(mask)+1;
	
	for(r=0; r < maxR; r++) {
		
		rmask = CreateImage(mask->ncols,mask->nrows);

		for(i=0;i<(mask->ncols*mask->nrows);i++) {
			if(mask->val[i] == r)
				rmask->val[i] = 1;
		}
		  
		  
		  sprintf(outputfile, "%s%d.fv", argv[3], r);

		  tic = Tic();
		  hist = BIC(cimg, rmask);
		  toc = Toc();
		  printf("\rCreateBIC for region %d in %f milliseconds\n",r, CTime(tic, toc));
		  WriteFileHistogram(hist,outputfile);
		  
		  DestroyImage(&rmask);
		  DestroyHistogram(&hist);
	}

  DestroyCImage(&cimg);
  DestroyImage(&mask);

  return(0);
}
