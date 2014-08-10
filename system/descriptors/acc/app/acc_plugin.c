#include "../include/libcolordescriptors.h"

void Extraction(char *img_path, char *fv_path)
{
  CImage *cimg=NULL;
  Histogram *acc=NULL;

  cimg = ReadCImage(img_path);

  acc = ACC(cimg);

  WriteFileHistogram(acc,fv_path);
  DestroyHistogram(&acc);
  DestroyCImage(&cimg);

}

void* LoadFV(char* fv_path) {
    return (void*) ReadFileHistogram(fv_path);
}

double Distance(void *fv1, void *fv2)
{
  Histogram *acc1=NULL;
  Histogram *acc2=NULL;
  double distance;

  acc1 = (Histogram*) fv1;
  acc2 = (Histogram*) fv2;

  distance = (double) L1Distance(acc1, acc2);

  //DestroyHistogram(&acc1);
  //DestroyHistogram(&acc2);

  return distance;
}
