#include "image.h"
#include "compressed.h"
#include "bzlib.h"

void WriteCompressedScene(Scene *scn, char *filename) {
  
  FILE   *f;
  //FILE   *in;
  BZFILE *b;
  //int     nBuf = 1024;
  //char    buf[1024];
  int    Imax, n, i;
  int    bzerror;
  char   *data;
  uchar  *data8;
  ushort *data16;

  //WriteScene(scn, "tmp.scn");
  //in = fopen("tmp.scn", "r");
  f = fopen( filename, "wb" );
  
  //if ( ( !f ) || ( !in ) ) {
  if ( !f ) {
    printf("Error on opening the compressed file %s\n", filename);
    return;
  }
  
  b = BZ2_bzWriteOpen( &bzerror, f, 9, 0, 30 );
  if (bzerror != BZ_OK) {
    BZ2_bzWriteClose ( &bzerror, b, 0, 0, 0 );
    fclose(f);
    printf("Error on opening the file %s\n", filename);
  }
  
  //while ( (bzerror == BZ_OK ) && (nBuf > 0) ) {
  if (bzerror == BZ_OK) {
    /* get data to write into buf, and set nBuf appropriately */
    data = (char*) calloc(200, sizeof(int));
    n = scn->xsize*scn->ysize*scn->zsize;
    sprintf(data,"SCN\n");
    sprintf(data,"%s%d %d %d\n",data,scn->xsize,scn->ysize,scn->zsize);
    sprintf(data,"%s%f %f %f\n",data,scn->dx,scn->dy,scn->dz);
    Imax = MaximumValue3(scn);
    if (Imax < 256) {
      printf("8bits\n");
      sprintf(data,"%s%d\n",data,8);
      BZ2_bzWrite ( &bzerror, b, data, strlen( data ) );
      data8 = AllocUCharArray(n);
      for (i=0; i < n; i++) 
	data8[i] = (uchar) scn->data[i];
      BZ2_bzWrite ( &bzerror, b, data8, n );
      free(data8);
    } else if (Imax < 65536) {
      printf("16bits\n");
      sprintf(data,"%s%d\n",data,16);
      BZ2_bzWrite ( &bzerror, b, data, strlen( data ) );
      data16 = AllocUShortArray(n);
      for (i=0; i < n; i++)
	data16[i] = (ushort) scn->data[i];
      BZ2_bzWrite ( &bzerror, b, data16, 2 * n );
      free(data16);
    } else {
      printf("32bits\n");
      sprintf(data,"%s%d\n",data,32);
      BZ2_bzWrite ( &bzerror, b, data, strlen( data ) );
      BZ2_bzWrite ( &bzerror, b, scn->data, 4 * n );
    }
    free(data);
    /*
      nBuf = fread(buf, sizeof(char), 1024, in);
      if (nBuf > 0)
      BZ2_bzWrite ( &bzerror, b, buf, nBuf );
      
      if (bzerror == BZ_IO_ERROR) { 
      break;
      }
    */
  }
  
  BZ2_bzWriteClose ( &bzerror, b, 0, 0, 0 );
  fclose(f);
  //fclose(in);

  if (bzerror == BZ_IO_ERROR) {
    printf("Error on writing to %s\n", filename);
  }
  
}


Scene *ReadCompressedScene(char *filename) {

  FILE * f, * out;
  BZFILE* b;
  int     nBuf;
  char    buf[1024];
  int     bzerror;
  Scene *scn = NULL;
  
  f = fopen ( filename, "r" );

  out = fopen("tmp.scn", "wb");

  if (( !f ) || ( !out )) {
    printf("Erro ao abrir arquivos!\n");
    return NULL;
  }
  
  b = BZ2_bzReadOpen ( &bzerror, f, 0, 0, NULL, 0 );

  if ( bzerror != BZ_OK ) {
    BZ2_bzReadClose ( &bzerror, b );
    fclose(f);
    fclose(out);
    printf("Erro ao abrir cena compactada!\n");
    return NULL;
  }
  
  bzerror = BZ_OK;
  nBuf = 1024;
  while ( bzerror == BZ_OK ) {
    nBuf = BZ2_bzRead ( &bzerror, b, buf, 1024);
    if ( ( bzerror == BZ_OK ) && (nBuf > 0) ) {
      fwrite(buf, sizeof(char), nBuf, out);
    }    
  }

  
  BZ2_bzReadClose ( &bzerror, b );
  fclose(f);
  fclose(out);
  
  scn = ReadScene("tmp.scn");

  return (scn);
}



// Read accordingly to the file extension.
Scene *ReadVolume(char *filename){
  Scene *scn=NULL;
  int s = strlen(filename);
  if(s>3 && strcasecmp(filename+s-3,"scn")==0)
    scn = ReadScene(filename);
  else if(s>7 && strcasecmp(filename+s-7,"scn.bz2")==0)
    scn = ReadCompressedScene(filename);
  else
    Warning("Unknown data type","ReadVolume");
  return scn;
}


// Write accordingly to the file extension.
void WriteVolume(Scene *scn, char *filename){
  int s = strlen(filename);
  if(s>3 && strcasecmp(filename+s-3,"scn")==0)
    WriteScene(scn, filename);
  else if(s>7 && strcasecmp(filename+s-7,"scn.bz2")==0)
    WriteCompressedScene(scn, filename);
  else
    Warning("Unknown data type","WriteVolume");
}



