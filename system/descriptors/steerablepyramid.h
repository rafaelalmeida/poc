#ifndef _STEERABLEPYRAMID_H_
#define _STEERABLEPYRAMID_H_

#include "dimage.h"
#include "image.h"
#include "spectrum.h"
#include "feature.h"
#include "common.h"
#include "featurevector.h"
#include <math.h>

/******************************************/
/**  Author: Thiago Vallin Spina         **/
/**  October, 2008                       **/
/**  Used in MO445                       **/
/**  thiago.spina@gmail.com              **/
/**                                      **/
/** Created according:                   **/
/** Combining global with local texture  **/
/** information for image retrieval      **/
/** applications - Montoya, J. et al.    **/
/******************************************/

//#define SPORIENTATIONS 6
//#define SPSCALES 4
#define SPORIENTATIONS 4
#define SPSCALES 2


/** Enum type used to define which filter will be created in
 *  CreateFilter
 **/
typedef enum filtertype{BandPass,LowPass} FilterType;


/** Filter creation and use **/

/** Extracts SP features from all features in feats.
 *  By default 4 scales and 6 orientations are used, but that can be
 *  changed by altering SPSCALES and SPORIENTATIONS.
 *
 *  @param feats Features* already extracted from the image (e.g. Lab, RGB, etc.).
 *
 *  @return a new Features* struct which will contain feats->nfeats*SCALES*SPORIENTATIONS features,
 *      because each existent feature in feats will be filtered using SPSCALES and SPORIENTATIONS
 **/
 
Features* SteerablePyramidFeats(Features* feats);
void WriteFileFeatureVector1D_bin(FeatureVector1D *fv,char *filename); //gera vetor em formato binario
FeatureVector1D *ReadFileFeatureVector1D_bin(char *filename);

#endif // _STEERABLEPYRAMID_H_
