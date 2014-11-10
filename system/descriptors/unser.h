#ifndef _UNSER_H_
#define _UNSER_H_

#include "image.h"
#include "histogram.h"

Histogram *Unser(Image *img, Image *msk);
int UnserDimensions();

#endif
