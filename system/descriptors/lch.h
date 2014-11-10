#ifndef _LCH_H_
#define _LCH_H_

#include "cimage.h"
#include "histogram.h"

Histogram *LCH(CImage *cimg, Image *mask);
int LCHDimensions();

#endif
