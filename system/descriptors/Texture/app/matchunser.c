#include "texturedescriptorslib.h"

int main(int argc, char** argv)
{
  timer *tic, *toc;
  Histogram *unser1=NULL;
  Histogram *unser2=NULL;
  ulong distance;

  if (argc != 3) {
    fprintf(stderr,"usage: matchgch <curve1> <curve2>\n");
    exit(-1);
  }  

  unser1 = ReadFileHistogram(argv[1]);
  unser2 = ReadFileHistogram(argv[2]);

  tic = Tic();
  distance=L1Distance(unser1, unser2);
  toc = Toc();
  printf("%ld\n",distance);

  DestroyHistogram(&unser1);
  DestroyHistogram(&unser2);

  return(0);
}
