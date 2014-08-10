#include "jpeg.h"

#include <jpeglib.h>

CImage *ReadJPEGFile(char *filename)
{
  CImage *cimg;
  FILE *fp;

  fp = fopen (filename, "rb");
  if (fp == NULL){
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(-1);
  }

  cimg = ReadJPEG(fp);
  fclose(fp);

  return(cimg);
}

CImage *ReadJPEG(FILE *fp)
{
  CImage *cimg;
  struct jpeg_decompress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPARRAY buffer;
  int row_stride;
  int p, c;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_decompress(&cinfo);

  jpeg_stdio_src(&cinfo, fp);

  (void) jpeg_read_header(&cinfo, TRUE);
  (void) jpeg_start_decompress(&cinfo);

  row_stride = cinfo.output_width * cinfo.output_components;
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

  cimg = CreateCImage(cinfo.output_height, cinfo.output_width);
  p = 0;
  while (cinfo.output_scanline < cinfo.output_height) {
    (void) jpeg_read_scanlines(&cinfo, buffer, 1);
    for (c = 0; c < cinfo.output_width; c++) {
      cimg->C[0]->val[p] = (int) buffer[0][c * cinfo.output_components + 0];
      cimg->C[1]->val[p] = (int) buffer[0][c * cinfo.output_components + 1];
      cimg->C[2]->val[p] = (int) buffer[0][c * cinfo.output_components + 2];
      p++;
    }
  }
  (void) jpeg_finish_decompress(&cinfo);
  jpeg_destroy_decompress(&cinfo);

  return(cimg);
}

void WriteJPEGFile(char *filename, CImage *cimg)
{
  FILE *fp;

  fp = fopen (filename, "wb");
  if (fp == NULL){
    fprintf(stderr,"Cannot open %s\n",filename);
    exit(-1);
  }

  WriteJPEG(fp, cimg);
  fclose(fp);
}

void WriteJPEG(FILE *fp, CImage *cimg)
{
  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;
  JSAMPARRAY buffer;
  int row_stride;
  int p, c;

  cinfo.err = jpeg_std_error(&jerr);
  jpeg_create_compress(&cinfo);

  jpeg_stdio_dest(&cinfo, fp);

  cinfo.image_width = cimg->C[0]->ncols;
  cinfo.image_height = cimg->C[0]->nrows;
  cinfo.input_components = 3;
  cinfo.in_color_space = JCS_RGB;

  jpeg_set_defaults(&cinfo);

  jpeg_start_compress(&cinfo, TRUE);

  row_stride = cinfo.image_width * cinfo.input_components;
  buffer = (*cinfo.mem->alloc_sarray)
		((j_common_ptr) &cinfo, JPOOL_IMAGE, row_stride, 1);

  p = 0;
  while (cinfo.next_scanline < cinfo.image_height) {
    for (c = 0; c < cinfo.image_width; c++) {
      buffer[0][c * cinfo.input_components + 0] = (uchar)cimg->C[0]->val[p];
      buffer[0][c * cinfo.input_components + 1] = (uchar)cimg->C[1]->val[p];
      buffer[0][c * cinfo.input_components + 2] = (uchar)cimg->C[2]->val[p];
      p++;
    }
    (void) jpeg_write_scanlines(&cinfo, buffer, 1);
  }
  jpeg_finish_compress(&cinfo);
  jpeg_destroy_compress(&cinfo);
}

