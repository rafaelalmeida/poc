#ifndef _GEOMETRIC2_H_
#define _GEOMETRIC2_H_

#include "scene.h"
#include "geometry.h"
#include "algebra.h"
#include "geometric3.h"

void transformacao2(Image *img, float IM[4][4], Image *img_transf);
void transformacaoDireta2(Image *img, float M[4][4], Image *img_transf);
void Rotacao2(float th, float T[4][4],Image *img ,Image *out);
void Translacao2(float dx,float dy,float T[4][4],Image *img, Image *out);
void RotTrans2(float th,float dx,float dy,float T[4][4],Image *img,Image *out);
#endif
