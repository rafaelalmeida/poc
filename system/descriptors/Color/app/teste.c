#include "colordescriptorslib.h"

// Cria uma imagem com os pixels da regiao
CImage* createRegionImage(int r, Image *mask, CImage* original) {
	
	int x, y, max_x, min_x, max_y, min_y, index, h, w, np, i;
	double meanR, meanG, meanB;
	CImage *cimg=NULL;
	
	min_x = mask->ncols-1;
	min_y = mask->nrows-1;
	max_x = 0;
	max_y = 0;
	
	meanR = 0;
	meanG = 0;
	meanB = 0;
	
	np=0;
	
	// Computing the bounding box
	for(y=0; y < mask->nrows; y++) {
			for(x=0; x < mask->ncols; x++) {
				
				index = y*(mask->ncols)+x;
				
				if(mask->val[index] == r) {
					
					if(x < min_x)
						min_x = x;
					if(y < min_y)
						min_y = y;
					if(x > max_x)
						max_x = x;
					if(y > max_y)
						max_y = y;
						
					meanR += (original->C[0]->val[index]);
					meanG += (original->C[1]->val[index]);
					meanB += (original->C[2]->val[index]);
					
					np++;
				}
			
			}
	}
	
	h = (max_y-min_y)+1;
	w = (max_x-min_x)+1;
	
	meanR /= np; 
	meanG /= np;
	meanB /= np;
	
	printf("Region %d \t x1=%d y1=%d x2=%d y2=%d", r, min_x, min_y, max_x, max_y);
	printf(" \t W=%d H=%d", w, h);
	printf("\tMean (R,G,B): (%d, %d, %d)\n", (int)meanR, (int)meanG, (int)meanB);
	
	// Construindo a nova imagem	
	cimg = CreateCImage(w,h);
	i=0;
	
	for(y=min_y; y <= max_y; y++) {
			for(x=min_x; x <= max_x; x++) {
				
				index = y*(mask->ncols)+x;
				
				if(mask->val[index] == r) {
					cimg->C[0]->val[i] = original->C[0]->val[index];
					cimg->C[1]->val[i] = original->C[1]->val[index];
					cimg->C[2]->val[i] = original->C[2]->val[index];
				} else {
//					cimg->C[0]->val[i] = meanR;
//					cimg->C[1]->val[i] = meanG;
//					cimg->C[2]->val[i] = meanB;
					cimg->C[0]->val[i] = 0;
					cimg->C[1]->val[i] = 0;
					cimg->C[2]->val[i] = 0;
				}
				
				i++;				
			}
	}
	
	return cimg;
}

int main(int argc, char** argv)
{
  timer *tic, *toc;
  CImage *cimg=NULL;
  Image *mask=NULL;
  CImage *rmask=NULL;
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
	
		rmask = createRegionImage(r, mask, cimg);
		
		sprintf(outputfile, "%s%d.ppm", argv[3], r);
		WriteCImage(rmask,outputfile);
		
/*		rmask = CreateImage(mask->ncols,mask->nrows);

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
		  
		  
		  DestroyHistogram(&hist);*/
		  
		  DestroyCImage(&rmask);
	}

  DestroyCImage(&cimg);
  DestroyImage(&mask);

  return(0);
}
