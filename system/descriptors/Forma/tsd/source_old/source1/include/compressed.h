
#ifndef _COMPRESSED_H_
#define _COMPRESSED_H_

#include "ift.h"

#define COMPRESSED_TYPE 0
#define NORMAL_TYPE 1

void   WriteCompressedScene(Scene *scn, 
			    char *filename);
Scene *ReadCompressedScene(char *filename);


// Read/Write accordingly to the file extension.
Scene *ReadVolume(char *filename);
void   WriteVolume(Scene *scn, char *filename);

#endif

