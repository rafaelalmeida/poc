#include "libtexturedescriptors.h"

double Distance(char *fv1_path, char *fv2_path)
{
  
  FeatureVector1D *steerable1=NULL;
  FeatureVector1D *steerable2=NULL;
  double distance;

  steerable1 = ReadFileFeatureVector1D_bin(fv1_path);
  steerable2 = ReadFileFeatureVector1D_bin(fv2_path);

  //distance = EuclideanDistance(steerable1, steerable2);
  distance = ManhattanDistance(steerable1, steerable2);

  DestroyFeatureVector1D(&steerable1);
  DestroyFeatureVector1D(&steerable2);

  return (distance);
}

int main(int argc, char** argv)
{
  double distance;

  if (argc != 3) {
    fprintf(stderr,"usage: sid_distance <fv1_path> <fv2_path>\n");
    exit(-1);
  }

  distance=Distance(argv[1], argv[2]);

  printf("%f\n",distance);

  return(0);
}
