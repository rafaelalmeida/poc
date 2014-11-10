#ifndef _BIC_H_
#define _BIC_H_

#include "cimage.h"
#include "histogram.h"

int BICDimensions();
Histogram *BIC(CImage *cimg, Image *mask);

#endif
