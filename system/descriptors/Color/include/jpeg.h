#ifndef _JPEG_H_
#define _JPEG_H_

#include "cimage.h"

CImage *ReadJPEGFile(char *filename);
CImage *ReadJPEG(FILE *fp);
void WriteJPEGFile(char *filename, CImage *cimg);
void WriteJPEG(FILE *fp, CImage *cimg);

#endif
