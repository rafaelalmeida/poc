#include "common.h"

char *AllocCharArray(int n)
{
  char *v=NULL;
  v = (char *) calloc(n,sizeof(char));
  if (v == NULL)
    Error(MSG1,"AllocCharArray");
  return(v);
}

uchar *AllocUCharArray(int n)
{
  uchar *v=NULL;
  v = (uchar *) calloc(n,sizeof(uchar));
  if (v == NULL)
    Error(MSG1,"AllocUCharArray");
  return(v);
}

ushort *AllocUShortArray(int n)
{
  ushort *v=NULL;
  v = (ushort *) calloc(n,sizeof(ushort));
  if (v == NULL)
    Error(MSG1,"AllocUShortArray");
  return(v);
}

ulong *AllocULongArray(int n)
{
  ulong *v=NULL;
  v = (ulong *) calloc(n,sizeof(ulong));
  if (v == NULL)
    Error(MSG1,"AllocULongArray");
  return(v);
}

bool *AllocBoolArray(int n)
{
  bool *v=NULL;
  v = (bool *) calloc(n,sizeof(bool));
  if (v == NULL)
    Error(MSG1,"AllocBoolArray");
  return(v);
}

int *AllocIntArray(int n)
{
  int *v=NULL;
  v = (int *) calloc(n,sizeof(int));
  if (v == NULL)
    Error(MSG1,"AllocIntArray");
  return(v);
}

float *AllocFloatArray(int n)
{
  float *v=NULL;
  v = (float *) calloc(n,sizeof(float));
  if (v == NULL)
    Error(MSG1,"AllocFloatArray");
  return(v);
}

double *AllocDoubleArray(int n)
{
  double *v=NULL;
  v = (double *) calloc(n,sizeof(double));
  if (v == NULL)
    Error(MSG1,"AllocDoubleArray");
  return(v);
}

void Error(char *msg,char *func){ /* It prints error message and exits
                                    the program. */
  fprintf(stderr,"Error:%s in %s\n",msg,func);
  exit(-1);
}

void Warning(char *msg,char *func){ /* It prints warning message and
                                       leaves the routine. */
 fprintf(stdout,"Warning:%s in %s\n",msg,func);

}

void Change(size_t w, void *a, void *b){/* It changes content between a and b */
  char tmp;
  char *pa, *pb;

  pa = (char *)a;
  pb = (char *)b;
  while (w--) {
    tmp = *pa;
    *pa++ = *pb;
    *pb++ = tmp;
  }
}

int NCFgets(char *s, int m, FILE *f) {
  while(fgets(s,m,f)!=NULL)
    if (s[0]!='#') return 1;
  return 0;
}

int RandomInteger (int low, int high){
  int k;
  double d;

  d = (double) random () / ((double) RAND_MAX);
  k = d * (high - low);
  return low + k;
}

int SafeMod(int a, int n)
{
	int r = a % n;

	return (r >= 0) ? r : n+r;
}

uchar ComputeNorm(float value)
{
  return((uchar)(255. * value));
}

uchar ComputeLog(float value)
{
  uchar result;
  
  value = 255. * value;
  if(value==0.)       result=0;
  else if(value<1.)   result=1;
  else if(value<2.)   result=2;
  else if(value<4.)   result=3;
  else if(value<8.)   result=4;
  else if(value<16.)  result=5;
  else if(value<32.)  result=6;
  else if(value<64.)  result=7;
  else if(value<128.) result=8;
  else                result=9;
  
  return(result);
}

Property *AllocPropertyArray(int n)
{
  Property *v=NULL;
  v = (Property *) calloc(n,sizeof(Property));
  if (v==NULL)
    Error(MSG1,"AllocPropertyArray");
  return(v);
}

double L2DoubleDistance(double *v1, double *v2, int size) {
    int i;
    double d=0.0;

    for (i=0; i<size; i++) {
        d += pow((v1[i]-v2[i]), 2);
    }
    return sqrt(d);
}

double *ReadFileDouble(char *filename, int size) {
    FILE *fp;
    double *fv=NULL;

    fv = AllocDoubleArray(size);

    if ((fp = fopen(filename, "rb")) == NULL) {
        printf("ERRO CRIANDO ARQUIVO DE FV\n");
        exit(0);
    }
    fread(fv, sizeof(double), size, fp);
    fclose(fp);

    return fv;
}

float *ReadFileFloat(char *filename, int size) {
    FILE *fp;
    float *fv=NULL;

    fv = AllocFloatArray(size);

    if ((fp = fopen(filename, "rb")) == NULL) {
        printf("ERRO CRIANDO ARQUIVO DE FV\n");
        exit(0);
    }
    fread(fv, sizeof(float), size, fp);
    fclose(fp);

    return fv;
}

void WriteFileDouble(char *filename, double *vet, int size) {
    FILE *fout;

    if ((fout = fopen(filename, "wb")) == NULL) {
        printf("ERRO CRIANDO ARQUIVO DE FV\n");
        exit(0);
    }
    fwrite(vet, sizeof(double), size, fout);
    fclose(fout);
}

void WriteFileFloat(char *filename, float *vet, int size) {
    FILE *fout;

    if ((fout = fopen(filename, "wb")) == NULL) {
        printf("ERRO CRIANDO ARQUIVO DE FV\n");
        exit(0);
    }
    fwrite(vet, sizeof(float), size, fout);
    fclose(fout);
}

double L1FloatDistance(float *v1, float *v2, int size) {
    int i;
    double d=0.0;

    for (i=0; i<size; i++) {
        d += (double)fabsf(v1[i]-v2[i]);
    }
    return d;
}

//Encontra o minimo valor entre 3
float min(float x, float y, float z) {

    if ( (x<=y) && (x<=z) ) {
        return x;
    } else if ( (y<=x) && (y<=z) ) {
        return y;
    } else if ( (z<=x) && (z<=y) ) {
        return z;
    }
    return -1;
}

//Encontra o maximo valor entre 3
float max(float x, float y, float z) {

    if ( (x>=y) && (x>=z) ) {
        return x;
    } else if ( (y>=x) && (y>=z) ) {
        return y;
    } else if ( (z>=x) && (z>=y) ) {
        return z;
    }
    return -1;
}
