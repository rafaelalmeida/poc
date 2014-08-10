LIB=./lib
INCLUDE=./include
BIN=./bin
SRC=./src
OBJ=./obj
DEMO=./demo
# FLAGS= -g -Wall -D _DEBUG
FLAGS=  -O3  -Wall -fPIC

libift: $(LIB)/libift.a
	echo "libift.a built..."

$(LIB)/libift.a: \
$(OBJ)/common.o \
$(OBJ)/comptime.o \
$(OBJ)/border.o \
$(OBJ)/heap.o \
$(OBJ)/gheap.o \
$(OBJ)/realheap.o \
$(OBJ)/queue.o \
$(OBJ)/sort.o \
$(OBJ)/set.o \
$(OBJ)/curve.o \
$(OBJ)/geometry.o \
$(OBJ)/algebra.o \
$(OBJ)/plane.o \
$(OBJ)/polynom.o \
$(OBJ)/image.o \
$(OBJ)/cimage.o \
$(OBJ)/dimage.o \
$(OBJ)/descriptor.o \
$(OBJ)/salience.o \
$(OBJ)/bas.o \
$(OBJ)/metrics.o \
$(OBJ)/scene.o \
$(OBJ)/spectrum.o \
$(OBJ)/adjacency.o \
$(OBJ)/adjacency3.o \
$(OBJ)/annimg.o \
$(OBJ)/annscn.o \
$(OBJ)/frame.o \
$(OBJ)/genift.o \
$(OBJ)/genift3.o \
$(OBJ)/matrix.o \
$(OBJ)/realmatrix.o \
$(OBJ)/mathematics.o \
$(OBJ)/mathematics3.o \
$(OBJ)/radiometric.o \
$(OBJ)/radiometric3.o \
$(OBJ)/filtering.o \
$(OBJ)/filtering3.o \
$(OBJ)/morphology.o \
$(OBJ)/morphology3.o \
$(OBJ)/segmentation.o \
$(OBJ)/segmentation3.o \
$(OBJ)/analysis.o \
$(OBJ)/analysis3.o \
$(OBJ)/context.o \
$(OBJ)/graphics.o \
$(OBJ)/shading.o \
$(OBJ)/shell.o \
$(OBJ)/normal.o \
$(OBJ)/raycasting.o \
$(OBJ)/splatting.o \
$(OBJ)/reslicing.o \
$(OBJ)/classify.o \
$(OBJ)/spline.o \
$(OBJ)/stackandfifo.o \
$(OBJ)/color.o \
$(OBJ)/cimg.o \
$(OBJ)/tensorscale.o \
$(OBJ)/complgraph.o \
$(OBJ)/prscene.o \
$(OBJ)/gqueue.o \
$(OBJ)/ctree.o \
$(OBJ)/mclouds.o \
$(OBJ)/interpolation.o \
$(OBJ)/interpolation3.o \
$(OBJ)/matching.o \
$(OBJ)/geometric3.o \
$(OBJ)/geometric2.o \
$(OBJ)/registration.o \
$(OBJ)/reslice3.o \
$(OBJ)/glist.o \
$(OBJ)/opf.o \
$(OBJ)/tira.o 
	ar csr $(LIB)/libift.a \
$(OBJ)/common.o \
$(OBJ)/comptime.o \
$(OBJ)/border.o \
$(OBJ)/heap.o \
$(OBJ)/gheap.o \
$(OBJ)/realheap.o \
$(OBJ)/queue.o \
$(OBJ)/sort.o \
$(OBJ)/set.o \
$(OBJ)/curve.o \
$(OBJ)/geometry.o \
$(OBJ)/algebra.o \
$(OBJ)/plane.o \
$(OBJ)/polynom.o \
$(OBJ)/image.o \
$(OBJ)/cimage.o \
$(OBJ)/dimage.o \
$(OBJ)/descriptor.o \
$(OBJ)/salience.o \
$(OBJ)/bas.o \
$(OBJ)/metrics.o \
$(OBJ)/scene.o \
$(OBJ)/spectrum.o \
$(OBJ)/adjacency.o \
$(OBJ)/adjacency3.o \
$(OBJ)/annimg.o \
$(OBJ)/annscn.o \
$(OBJ)/frame.o \
$(OBJ)/genift.o \
$(OBJ)/genift3.o \
$(OBJ)/matrix.o \
$(OBJ)/realmatrix.o \
$(OBJ)/mathematics.o \
$(OBJ)/mathematics3.o \
$(OBJ)/radiometric.o \
$(OBJ)/radiometric3.o \
$(OBJ)/filtering.o \
$(OBJ)/filtering3.o \
$(OBJ)/morphology.o \
$(OBJ)/morphology3.o \
$(OBJ)/segmentation.o \
$(OBJ)/segmentation3.o \
$(OBJ)/analysis.o \
$(OBJ)/analysis3.o \
$(OBJ)/context.o \
$(OBJ)/graphics.o \
$(OBJ)/shading.o \
$(OBJ)/shell.o \
$(OBJ)/normal.o \
$(OBJ)/raycasting.o \
$(OBJ)/splatting.o \
$(OBJ)/reslicing.o \
$(OBJ)/classify.o \
$(OBJ)/spline.o \
$(OBJ)/stackandfifo.o \
$(OBJ)/color.o \
$(OBJ)/cimg.o \
$(OBJ)/tensorscale.o \
$(OBJ)/complgraph.o \
$(OBJ)/prscene.o \
$(OBJ)/gqueue.o \
$(OBJ)/ctree.o \
$(OBJ)/mclouds.o \
$(OBJ)/interpolation.o \
$(OBJ)/interpolation3.o \
$(OBJ)/matching.o \
$(OBJ)/geometric3.o \
$(OBJ)/geometric2.o \
$(OBJ)/registration.o \
$(OBJ)/reslice3.o \
$(OBJ)/glist.o \
$(OBJ)/opf.o \
$(OBJ)/tira.o 

$(OBJ)/glist.o: $(SRC)/glist.c
	gcc $(FLAGS) -c $(SRC)/glist.c -I$(INCLUDE) \
	-o $(OBJ)/glist.o 

$(OBJ)/common.o: $(SRC)/common.c
	gcc $(FLAGS) -c $(SRC)/common.c -I$(INCLUDE) \
	-o $(OBJ)/common.o 

$(OBJ)/comptime.o: $(SRC)/comptime.c
	gcc $(FLAGS) -c $(SRC)/comptime.c -I$(INCLUDE) \
	-o $(OBJ)/comptime.o 

$(OBJ)/border.o: $(SRC)/border.c
	gcc $(FLAGS) -c $(SRC)/border.c -I$(INCLUDE) \
	-o $(OBJ)/border.o 

$(OBJ)/heap.o: $(SRC)/heap.c
	gcc $(FLAGS) -c $(SRC)/heap.c -I$(INCLUDE) \
	-o $(OBJ)/heap.o 

$(OBJ)/gheap.o: $(SRC)/gheap.c
	gcc $(FLAGS) -c $(SRC)/gheap.c -I$(INCLUDE) \
	-o $(OBJ)/gheap.o 

$(OBJ)/realheap.o: $(SRC)/realheap.c
	gcc $(FLAGS) -c $(SRC)/realheap.c -I$(INCLUDE) \
	-o $(OBJ)/realheap.o

$(OBJ)/queue.o: $(SRC)/queue.c
	gcc $(FLAGS) -c $(SRC)/queue.c -I$(INCLUDE) \
	-o $(OBJ)/queue.o 

$(OBJ)/sort.o: $(SRC)/sort.c
	gcc $(FLAGS) -c $(SRC)/sort.c -I$(INCLUDE) \
	-o $(OBJ)/sort.o 

$(OBJ)/set.o: $(SRC)/set.c
	gcc $(FLAGS) -c $(SRC)/set.c -I$(INCLUDE) \
	-o $(OBJ)/set.o 

$(OBJ)/curve.o: $(SRC)/curve.c
	gcc $(FLAGS) -c $(SRC)/curve.c -I$(INCLUDE) \
	-o $(OBJ)/curve.o 

$(OBJ)/geometry.o: $(SRC)/geometry.c
	gcc $(FLAGS) -c $(SRC)/geometry.c -I$(INCLUDE) \
	-o $(OBJ)/geometry.o 

$(OBJ)/algebra.o: $(SRC)/algebra.c
	gcc $(FLAGS) -c $(SRC)/algebra.c -I$(INCLUDE) \
	-o $(OBJ)/algebra.o 

$(OBJ)/plane.o: $(SRC)/plane.c
	gcc $(FLAGS) -c $(SRC)/plane.c -I$(INCLUDE) \
	-o $(OBJ)/plane.o 

$(OBJ)/polynom.o: $(SRC)/polynom.c
	gcc $(FLAGS) -c $(SRC)/polynom.c -I$(INCLUDE) \
	-o $(OBJ)/polynom.o 

$(OBJ)/image.o: $(SRC)/image.c
	gcc $(FLAGS) -c $(SRC)/image.c -I$(INCLUDE) \
	-o $(OBJ)/image.o 

$(OBJ)/cimage.o: $(SRC)/cimage.c
	gcc $(FLAGS) -c $(SRC)/cimage.c -I$(INCLUDE) \
	-o $(OBJ)/cimage.o 

$(OBJ)/dimage.o: $(SRC)/dimage.c
	gcc $(FLAGS) -c $(SRC)/dimage.c -I$(INCLUDE) \
	-o $(OBJ)/dimage.o 

$(OBJ)/scene.o: $(SRC)/scene.c
	gcc $(FLAGS) -c $(SRC)/scene.c -I$(INCLUDE) \
	-o $(OBJ)/scene.o 

$(OBJ)/spectrum.o: $(SRC)/spectrum.c
	gcc $(FLAGS) -c $(SRC)/spectrum.c -I$(INCLUDE) \
	-o $(OBJ)/spectrum.o 

$(OBJ)/adjacency.o: $(SRC)/adjacency.c
	gcc $(FLAGS) -c $(SRC)/adjacency.c -I$(INCLUDE) \
	-o $(OBJ)/adjacency.o 

$(OBJ)/adjacency3.o: $(SRC)/adjacency3.c
	gcc $(FLAGS) -c $(SRC)/adjacency3.c -I$(INCLUDE) \
	-o $(OBJ)/adjacency3.o 

$(OBJ)/annimg.o: $(SRC)/annimg.c
	gcc $(FLAGS) -c $(SRC)/annimg.c -I$(INCLUDE) \
	-o $(OBJ)/annimg.o 

$(OBJ)/annscn.o: $(SRC)/annscn.c
	gcc $(FLAGS) -c $(SRC)/annscn.c -I$(INCLUDE) \
	-o $(OBJ)/annscn.o 

$(OBJ)/frame.o: $(SRC)/frame.c
	gcc $(FLAGS) -c $(SRC)/frame.c -I$(INCLUDE) \
	-o $(OBJ)/frame.o 

$(OBJ)/mathematics.o: $(SRC)/mathematics.c
	gcc $(FLAGS) -c $(SRC)/mathematics.c -I$(INCLUDE) \
	-o $(OBJ)/mathematics.o 

$(OBJ)/mathematics3.o: $(SRC)/mathematics3.c
	gcc $(FLAGS) -c $(SRC)/mathematics3.c -I$(INCLUDE) \
	-o $(OBJ)/mathematics3.o 

$(OBJ)/radiometric.o: $(SRC)/radiometric.c
	gcc $(FLAGS) -c $(SRC)/radiometric.c -I$(INCLUDE) \
	-o $(OBJ)/radiometric.o 

$(OBJ)/radiometric3.o: $(SRC)/radiometric3.c
	gcc $(FLAGS) -c $(SRC)/radiometric3.c -I$(INCLUDE) \
	-o $(OBJ)/radiometric3.o 

$(OBJ)/morphology.o: $(SRC)/morphology.c
	gcc $(FLAGS) -c $(SRC)/morphology.c -I$(INCLUDE) \
	-o $(OBJ)/morphology.o 

$(OBJ)/morphology3.o: $(SRC)/morphology3.c
	gcc $(FLAGS) -c $(SRC)/morphology3.c -I$(INCLUDE) \
	-o $(OBJ)/morphology3.o 

$(OBJ)/segmentation.o: $(SRC)/segmentation.c
	gcc $(FLAGS) -c $(SRC)/segmentation.c -I$(INCLUDE) \
	-o $(OBJ)/segmentation.o 

$(OBJ)/segmentation3.o: $(SRC)/segmentation3.c
	gcc $(FLAGS) -c $(SRC)/segmentation3.c -I$(INCLUDE) \
	-o $(OBJ)/segmentation3.o 

$(OBJ)/analysis.o: $(SRC)/analysis.c
	gcc $(FLAGS) -c $(SRC)/analysis.c -I$(INCLUDE) \
	-o $(OBJ)/analysis.o 

$(OBJ)/analysis3.o: $(SRC)/analysis3.c
	gcc $(FLAGS) -c $(SRC)/analysis3.c -I$(INCLUDE) \
	-o $(OBJ)/analysis3.o 

$(OBJ)/filtering.o: $(SRC)/filtering.c
	gcc $(FLAGS) -c $(SRC)/filtering.c -I$(INCLUDE) \
	-o $(OBJ)/filtering.o 

$(OBJ)/filtering3.o: $(SRC)/filtering3.c
	gcc $(FLAGS) -c $(SRC)/filtering3.c -I$(INCLUDE) \
	-o $(OBJ)/filtering3.o 

$(OBJ)/context.o: $(SRC)/context.c
	gcc $(FLAGS) -c $(SRC)/context.c -I$(INCLUDE) \
	-o $(OBJ)/context.o 

$(OBJ)/graphics.o: $(SRC)/graphics.c
	gcc $(FLAGS) -c $(SRC)/graphics.c -I$(INCLUDE) \
	-o $(OBJ)/graphics.o 

$(OBJ)/shading.o: $(SRC)/shading.c
	gcc $(FLAGS) -c $(SRC)/shading.c -I$(INCLUDE) \
	-o $(OBJ)/shading.o 

$(OBJ)/shell.o: $(SRC)/shell.c
	gcc $(FLAGS) -c $(SRC)/shell.c -I$(INCLUDE) \
	-o $(OBJ)/shell.o 

$(OBJ)/normal.o: $(SRC)/normal.c
	gcc $(FLAGS) -c $(SRC)/normal.c -I$(INCLUDE) \
	-o $(OBJ)/normal.o 

$(OBJ)/raycasting.o: $(SRC)/raycasting.c
	gcc $(FLAGS) -c $(SRC)/raycasting.c -I$(INCLUDE) \
	-o $(OBJ)/raycasting.o 

$(OBJ)/splatting.o: $(SRC)/splatting.c
	gcc $(FLAGS) -c $(SRC)/splatting.c -I$(INCLUDE) \
	-o $(OBJ)/splatting.o 

$(OBJ)/reslicing.o: $(SRC)/reslicing.c
	gcc $(FLAGS) -c $(SRC)/reslicing.c -I$(INCLUDE) \
	-o $(OBJ)/reslicing.o 

$(OBJ)/genift.o: $(SRC)/genift.c
	gcc $(FLAGS) -c $(SRC)/genift.c -I$(INCLUDE) \
	-o $(OBJ)/genift.o 

$(OBJ)/genift3.o: $(SRC)/genift3.c
	gcc $(FLAGS) -c $(SRC)/genift3.c -I$(INCLUDE) \
	-o $(OBJ)/genift3.o 

$(OBJ)/descriptor.o: $(SRC)/descriptor.c
	gcc $(FLAGS) -c $(SRC)/descriptor.c -I$(INCLUDE) \
	-o $(OBJ)/descriptor.o 

$(OBJ)/salience.o: $(SRC)/salience.c
	gcc $(FLAGS) -c $(SRC)/salience.c -I$(INCLUDE) \
	-o $(OBJ)/salience.o 

$(OBJ)/bas.o: $(SRC)/bas.c
	gcc $(FLAGS) -c $(SRC)/bas.c -I$(INCLUDE) \
	-o $(OBJ)/bas.o

$(OBJ)/metrics.o: $(SRC)/metrics.c
	gcc $(FLAGS) -c $(SRC)/metrics.c -I$(INCLUDE) \
	-o $(OBJ)/metrics.o

$(OBJ)/matrix.o: $(SRC)/matrix.c
	gcc $(FLAGS) -c $(SRC)/matrix.c -I$(INCLUDE) \
	-o $(OBJ)/matrix.o 

$(OBJ)/realmatrix.o: $(SRC)/realmatrix.c
	gcc $(FLAGS) -c $(SRC)/realmatrix.c -I$(INCLUDE) \
	-o $(OBJ)/realmatrix.o 

$(OBJ)/classify.o: $(SRC)/classify.c
	gcc $(FLAGS) -c $(SRC)/classify.c -I$(INCLUDE) \
	-o $(OBJ)/classify.o

$(OBJ)/spline.o: $(SRC)/spline.c
	gcc $(FLAGS) -c $(SRC)/spline.c -I$(INCLUDE) \
	-o $(OBJ)/spline.o

$(OBJ)/stackandfifo.o: $(SRC)/stackandfifo.c
	gcc $(FLAGS) -c $(SRC)/stackandfifo.c -I$(INCLUDE) \
	-o $(OBJ)/stackandfifo.o

$(OBJ)/cimg.o: $(SRC)/cimg.c
	gcc $(FLAGS) -c $(SRC)/cimg.c -I$(INCLUDE) \
	-o $(OBJ)/cimg.o

$(OBJ)/color.o: $(SRC)/color.c
	gcc $(FLAGS) -c $(SRC)/color.c -I$(INCLUDE) \
	-o $(OBJ)/color.o

$(OBJ)/tensorscale.o: $(SRC)/tensorscale.c
	gcc $(FLAGS) -c $(SRC)/tensorscale.c -I$(INCLUDE) \
	-o $(OBJ)/tensorscale.o

$(OBJ)/complgraph.o: $(SRC)/complgraph.c
	gcc $(FLAGS) -c $(SRC)/complgraph.c -I$(INCLUDE) \
	-o $(OBJ)/complgraph.o

$(OBJ)/prscene.o: $(SRC)/prscene.c
	gcc $(FLAGS) -c $(SRC)/prscene.c -I$(INCLUDE) \
	-o $(OBJ)/prscene.o

$(OBJ)/gqueue.o: $(SRC)/gqueue.c
	gcc $(FLAGS) -c $(SRC)/gqueue.c -I$(INCLUDE) \
	-o $(OBJ)/gqueue.o

$(OBJ)/ctree.o: $(SRC)/ctree.c
	gcc $(FLAGS) -c $(SRC)/ctree.c -I$(INCLUDE) \
	-o $(OBJ)/ctree.o

$(OBJ)/mclouds.o: $(SRC)/mclouds.c
	gcc $(FLAGS) -c $(SRC)/mclouds.c -I$(INCLUDE) \
	-o $(OBJ)/mclouds.o

$(OBJ)/interpolation.o: $(SRC)/interpolation.c
	gcc $(FLAGS) -c $(SRC)/interpolation.c -I$(INCLUDE) \
	-o $(OBJ)/interpolation.o

$(OBJ)/interpolation3.o: $(SRC)/interpolation3.c
	gcc $(FLAGS) -c $(SRC)/interpolation3.c -I$(INCLUDE) \
	-o $(OBJ)/interpolation3.o
	
$(OBJ)/matching.o: $(SRC)/matching.c
	gcc $(FLAGS) -c $(SRC)/matching.c -I$(INCLUDE) \
	-o $(OBJ)/matching.o

$(OBJ)/geometric3.o: $(SRC)/geometric3.c
	gcc $(FLAGS) -c $(SRC)/geometric3.c -I$(INCLUDE) \
	-o $(OBJ)/geometric3.o
	
$(OBJ)/geometric2.o: $(SRC)/geometric2.c
	gcc $(FLAGS) -c $(SRC)/geometric2.c -I$(INCLUDE) \
	-o $(OBJ)/geometric2.o

$(OBJ)/registration.o: $(SRC)/registration.c
	gcc $(FLAGS) -c $(SRC)/registration.c -I$(INCLUDE) \
	-o $(OBJ)/registration.o

$(OBJ)/reslice3.o: $(SRC)/reslice3.c
	gcc $(FLAGS) -c $(SRC)/reslice3.c -I$(INCLUDE) \
	-o $(OBJ)/reslice3.o

$(OBJ)/opf.o: $(SRC)/opf.c
	gcc $(FLAGS) -c $(SRC)/opf.c -I$(INCLUDE) \
	-o $(OBJ)/opf.o

$(OBJ)/tira.o: $(SRC)/tira.c
	gcc $(FLAGS) -c $(SRC)/tira.c -I$(INCLUDE) \
	-o $(OBJ)/tira.o

BZIP2_DIR = ../bzip2
BZIP2FLAGS  = -I$(BZIP2_DIR)

$(OBJ)/compressed.o: $(SRC)/compressed.c
	gcc $(FLAGS) $(BZIP2FLAGS) -c $(SRC)/compressed.c \
	-I$(INCLUDE) -o $(OBJ)/compressed.o


compressed: libift \
$(OBJ)/compressed.o
	cp $(LIB)/libift.a $(LIB)/libiftc.a; ar sr $(LIB)/libiftc.a \
$(OBJ)/compressed.o


clean: 
	rm $(LIB)/lib*.a; rm $(OBJ)/*.o; rm $(BIN)/*	




