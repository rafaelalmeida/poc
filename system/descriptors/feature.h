#ifndef _FEATURES_H_
#define _FEATURES_H_

#include "dimage.h"
#include "spectrum.h"
#include "common.h"

typedef struct _felem {
  float *feat;
} FElem;

typedef struct _features {
  FElem *elem;
  int  nfeats;
  int  nelems;
  int  nrows,ncols;
  int  Imax;
} Features;

  typedef struct{
    Spectrum ** m_pSpectralFilterBank;
    Features * m_pFeatures;
    int orientations;
    int scales;
  }Filterbank;

/// Allocates memory for features
Features* CreateFeatures(int ncols, int nrows, int nfeats);

void DestroyFeatures(Features **f);

/** Gets the given feature for all the image pixels
 *  Returns a double image with the features
 **/
DImage* GetFeature(Features* feats, int index);

int SetFeature(Features* feats, int index, DImage* imagefeats);

/** Concatenates 2 features.
 * If either features are NULL, or they are of different sizes,
 * it will be returned NULL
 * NOTE: The features are not deallocated
 **/
Features* ConcatFeatures(Features* f1, Features* f2);

/** Used to return the image(Features*) to its original size,
 *  after it is power-of-2-sided by the Steerable
 *  Pyramid filter. Actually, the image is trimmed until
 *  it becomes ncols width and nrows height.
 *
 *  Created for Steerable Pyramid filtering.
 **/
Features* RemoveZeroes(Features* feats,int ncols,int nrows);

void FNormalizeFeatures(Features *f);

#endif
