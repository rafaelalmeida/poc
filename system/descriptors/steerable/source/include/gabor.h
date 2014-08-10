#ifndef _GABOR_H_
#define _GABOR_H_

#include <assert.h>
#include "dimage.h"
#include "image.h"
#include "spectrum.h"
#include "feature.h"
#include "common.h"

  #ifndef NSZE
    // number of scales during image decomposition.
    #define NSZE 128
  #endif

  #ifndef NSCALES
    // number of scales during image decomposition.
    #define NSCALES 4
  #endif

  #ifndef NORIENTATIONS
    // number of orientations during image decomposition.
    #define NORIENTATIONS 6
  #endif

  #ifndef DUH
    // highest frequency.
    #define DUH 0.4
  #endif

  #ifndef DUL
    // lowest frequency.
    #define DUL 0.05
  #endif
  
  typedef struct{
      DImage * m_pReal;
      DImage  * m_pImag;
      int m_nSze;
  }gabor;

  Features* GaborFeats(Features* feats);

#endif
