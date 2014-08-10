#include "colordescriptorslib.h"

int main(int argc, char** argv)
{
  timer *tic, *toc;
  Histogram *ccv1=NULL;
  Histogram *ccv2=NULL;
  ulong distance;

  if (argc != 3) {
    fprintf(stderr,"usage: matchccv <curve1> <curve2>\n");
    exit(-1);
  }  

  ccv1 = ReadFileHistogram(argv[1]);
  ccv2 = ReadFileHistogram(argv[2]);

  tic = Tic();
  distance=L1Distance(ccv1, ccv2);
  toc = Toc();
  printf("%ld\n",distance);

  DestroyHistogram(&ccv1);
  DestroyHistogram(&ccv2);

  return(0);
}
