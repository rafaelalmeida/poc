#include "colordescriptorslib.h"

int main(int argc, char** argv)
{
  timer *tic, *toc;
  Histogram *lch1=NULL;
  Histogram *lch2=NULL;
  ulong distance;

  if (argc != 3) {
    fprintf(stderr,"usage: matchlch <curve1> <curve2>\n");
    exit(-1);
  }  

  lch1 = ReadFileHistogram(argv[1]);
  lch2 = ReadFileHistogram(argv[2]);

  tic = Tic();
  distance=L1Distance(lch1, lch2);
  toc = Toc();
  printf("%ld\n",distance);

  DestroyHistogram(&lch1);
  DestroyHistogram(&lch2);

  return(0);
}
