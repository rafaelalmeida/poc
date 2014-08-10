#include "unser.h"

#include "common.h"
#include "adjacency.h"

void ComputeHistograms(Image *img, Image *msk,
                       float sum[4][511], float dif[4][511])
{
  ulong x, y, p, q;
  ulong npoints;
  AdjRel *A;
  int i, j;
  Pixel v;
  
  A = Circular(1.5);

  for (i=0; i<4; i++)
    for (j=0; j<=510; j++){
      sum[i][j] = 0.0;
      dif[i][j] = 0.0;
    }

  npoints = 0;
  for(y=1; y<img->nrows-1; y++){
    for(x=1; x<img->ncols-1; x++){
      p = x + img->tbrow[y];
      if (msk->val[p]) {
        for (i=1; i <= A->n>>1; i++){
          v.x = x + A->dx[i];
          v.y = y + A->dy[i];
          q = v.x + img->tbrow[v.y];
          if (ValidPixel(img,v.x,v.y) && msk->val[q]){
            sum[i - 1][img->val[p] + img->val[q]] += 1.0;
            dif[i - 1][img->val[p] - img->val[q] + 255] += 1.0;
	  }
        }
        npoints++;
      }
    }
  }
  DestroyAdjRel(&A);

  for (i=0; i<4; i++)
    for (j=0; j<=510; j++){
      sum[i][j] /= (float) npoints;
      dif[i][j] /= (float) npoints;
    }
}

float Mean(float s[511])
{
  float mean;
  int i;

  mean = 0.0;
  for (i=0; i<=510; i++)
    mean += i * s[i];
  mean *= 0.5;

  return(mean);
}

float Contrast(float d[511])
{
  float contrast;
  int j;

  contrast = 0.0;
  for (j=-255; j<=255; j++)
    contrast += j * j * d[j+255];

  return(contrast);
}

float Correlation(float s[511], float mean, float contrast)
{
  float correlation, aux;
  int i;

  correlation = 0.0;
  for (i=0; i<=510; i++){
    aux = i - 2.0 * mean;
    correlation += aux * aux * s[i];
  }
  correlation -= contrast;
  correlation *= 0.5;

  return(correlation);
}

float Energy(float s[511], float d[511])
{
  float energy_s, energy_d;
  int i;

  energy_s = 0.0;
  energy_d = 0.0;
  for (i=0; i<=510; i++){
    energy_s += s[i] * s[i];
    energy_d += d[i] * d[i];
  }

  return(energy_s * energy_d);
}

float Entropy(float s[511], float d[511])
{
  float entropy_s, entropy_d;
  int i;

  entropy_s = 0.0;
  entropy_d = 0.0;
  for (i=0; i<=510; i++){
    if (s[i] > 0.0)
      entropy_s += s[i] * log10(s[i]);
    if (d[i] > 0.0)
      entropy_d += d[i] * log10(d[i]);
  }

  return(- entropy_s - entropy_d);
}

float Homogeneity(float d[511])
{
  float homogeneity;
  int j;

  homogeneity = 0.0;
  for (j=-255; j<=255; j++)
    homogeneity += (1.0 / (1.0 + (float)(j * j))) * d[j+255];

  return(homogeneity);
}

float MaximalProbability(float s[511])
{
  float max;
  int i;

  max = 0.0;
  for (i=0; i<=510; i++)
    if (max < s[i])
      max = s[i];

  return(max);
}

float StandardDeviation(float contrast, float correlation)
{
  return(sqrt(correlation+contrast));
}

Histogram *Unser(Image *img, Image *msk)
{
  Histogram *h = NULL;
  float sum[4][511], dif[4][511];
  float mean, contrast, correlation;
  int i;
  
  ComputeHistograms(img, msk, sum, dif);

  h = CreateHistogram(32);
  for (i=0; i<4; i++){
    mean = Mean(sum[i]);
    h->v[i * 8 + 0] = (uchar)(mean);
    contrast = Contrast(dif[i]);
    h->v[i * 8 + 1] = (uchar)(contrast / 255.0);
    correlation = Correlation(sum[i], mean, contrast);
    h->v[i * 8 + 2] = (uchar)((correlation + 32512.5) / 255.0);
    h->v[i * 8 + 3] = (uchar)(255.0 * Energy(sum[i], dif[i]));
    h->v[i * 8 + 4] = (uchar)(255.0 * Entropy(sum[i], dif[i]) / 5.4168418);
    h->v[i * 8 + 5] = (uchar)(255.0 * Homogeneity(dif[i]));
    h->v[i * 8 + 6] = (uchar)(255.0 * MaximalProbability(sum[i]));
    h->v[i * 8 + 7] = (uchar)(sqrt(2)*StandardDeviation(contrast, correlation));
  }

  return(h);
}

