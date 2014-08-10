#include "colordescriptorslib.h"

int main(int argc, char** argv)
{
  timer *tic, *toc;
  Histogram *bic1=NULL;
  Histogram *bic2=NULL;
  ulong distance;

  if (argc != 3) {
    fprintf(stderr,"usage: matchbic <curve1> <curve2>\n");
    exit(-1);
  }  

  bic1 = ReadFileHistogram(argv[1]);
  bic2 = ReadFileHistogram(argv[2]);

  tic = Tic();
  distance=L1Distance(bic1, bic2);
  toc = Toc();
  printf("%ld\n",distance);

  DestroyHistogram(&bic1);
  DestroyHistogram(&bic2);

  return(0);
}
