#include "colordescriptorslib.h"

int main(int argc, char** argv)
{
  timer *tic, *toc;
  Histogram *gch1=NULL;
  Histogram *gch2=NULL;
  ulong distance;

  if (argc != 3) {
    fprintf(stderr,"usage: matchgch <curve1> <curve2>\n");
    exit(-1);
  }  

  gch1 = ReadFileHistogram(argv[1]);
  gch2 = ReadFileHistogram(argv[2]);

  tic = Tic();
  distance=L1Distance(gch1, gch2);
  toc = Toc();
  printf("%ld\n",distance);

  DestroyHistogram(&gch1);
  DestroyHistogram(&gch2);

  return(0);
}
