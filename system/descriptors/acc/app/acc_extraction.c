#include "libcolordescriptors.h"

int main(int argc, char** argv)
{
  timer *tic, *toc;
  CImage *cimg=NULL;
  Image *mask=NULL;
  Image *rmask=NULL;
  Histogram **acc=NULL;
  int i, r, maxR;
  char outputfile[100];
  
  if (argc != 4) {
    fprintf(stderr,"usage:  acc_extraction <image> <mask> <curve>\n");
    exit(-1);
  }  
  
	cimg = ReadCImage(argv[1]);
	mask = ReadImage(argv[2]);
  
	maxR = MaximumValue(mask)+1;
	
	acc = (Histogram**) calloc(maxR, sizeof(Histogram*));
	
	for(r=0; r < maxR; r++) {
		
		acc[r] = NULL;
		
		rmask = CreateImage(mask->ncols,mask->nrows);

		for(i=0;i<(mask->ncols*mask->nrows);i++) {
			if(mask->val[i] == r)
				rmask->val[i] = 1;
		}
		  
		  sprintf(outputfile, "%s%d.fv", argv[3], r);

		  tic = Tic();
		  acc[r] = ACC(cimg, rmask);
		  toc = Toc();
		  
		  printf("\rCreate ACC for region %d in %f milliseconds\n",r, CTime(tic, toc));
		  WriteFileHistogram(acc[r],outputfile);
		  
		  sprintf(outputfile, "%d%s.pgm", r, argv[3]);
		  //WriteImage(rmask,outputfile);		  
		  
		  DestroyImage(&rmask);
		  //DestroyHistogram(&acc[r]);
	}

  //free(acc);
  DestroyCImage(&cimg);
  DestroyImage(&mask);

  return(0);
}
