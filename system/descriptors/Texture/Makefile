LIB=./lib
INCLUDE=./include
SRC=./src
OBJ=./obj
FLAGS= -O2 -Wall

libtexturedescriptors: $(LIB)/libtexturedescriptors.a
	echo "libtexturedescriptors.a built..."

$(LIB)/libtexturedescriptors.a: \
$(OBJ)/common.o \
$(OBJ)/comptime.o \
$(OBJ)/image.o \
$(OBJ)/cimage.o \
$(OBJ)/adjacency.o \
$(OBJ)/mathematics.o \
$(OBJ)/color.o \
$(OBJ)/featurevector.o \
$(OBJ)/metrics.o \
$(OBJ)/histogram.o \
$(OBJ)/unser.o \
$(OBJ)/jpeg.o 

	ar csr $(LIB)/libtexturedescriptors.a \
$(OBJ)/common.o \
$(OBJ)/comptime.o \
$(OBJ)/image.o \
$(OBJ)/cimage.o \
$(OBJ)/adjacency.o \
$(OBJ)/mathematics.o \
$(OBJ)/color.o \
$(OBJ)/featurevector.o \
$(OBJ)/metrics.o \
$(OBJ)/histogram.o \
$(OBJ)/unser.o \
$(OBJ)/jpeg.o 

$(OBJ)/common.o: $(SRC)/common.c
	gcc $(FLAGS) -c $(SRC)/common.c -I$(INCLUDE) \
	-o $(OBJ)/common.o 

$(OBJ)/comptime.o: $(SRC)/comptime.c
	gcc $(FLAGS) -c $(SRC)/comptime.c -I$(INCLUDE) \
	-o $(OBJ)/comptime.o 

$(OBJ)/image.o: $(SRC)/image.c
	gcc $(FLAGS) -c $(SRC)/image.c -I$(INCLUDE) \
	-o $(OBJ)/image.o 

$(OBJ)/cimage.o: $(SRC)/cimage.c
	gcc $(FLAGS) -c $(SRC)/cimage.c -I$(INCLUDE) \
	-o $(OBJ)/cimage.o 

$(OBJ)/adjacency.o: $(SRC)/adjacency.c
	gcc $(FLAGS) -c $(SRC)/adjacency.c -I$(INCLUDE) \
	-o $(OBJ)/adjacency.o 

$(OBJ)/mathematics.o: $(SRC)/mathematics.c
	gcc $(FLAGS) -c $(SRC)/mathematics.c -I$(INCLUDE) \
	-o $(OBJ)/mathematics.o 

$(OBJ)/color.o: $(SRC)/color.c
	gcc $(FLAGS) -c $(SRC)/color.c -I$(INCLUDE) \
	-o $(OBJ)/color.o

$(OBJ)/featurevector.o: $(SRC)/featurevector.c
	gcc $(FLAGS) -c $(SRC)/featurevector.c -I$(INCLUDE) \
	-o $(OBJ)/featurevector.o

$(OBJ)/metrics.o: $(SRC)/metrics.c
	gcc $(FLAGS) -c $(SRC)/metrics.c -I$(INCLUDE) \
	-o $(OBJ)/metrics.o

$(OBJ)/histogram.o: $(SRC)/histogram.c
	gcc $(FLAGS) -c $(SRC)/histogram.c -I$(INCLUDE) \
	-o $(OBJ)/histogram.o

$(OBJ)/unser.o: $(SRC)/unser.c
	gcc $(FLAGS) -c $(SRC)/unser.c -I$(INCLUDE) \
	-o $(OBJ)/unser.o

$(OBJ)/jpeg.o: $(SRC)/jpeg.c
	gcc $(FLAGS) -c $(SRC)/jpeg.c -I$(INCLUDE) \
	-o $(OBJ)/jpeg.o

clean: 
	rm $(LIB)/lib*.a; rm $(OBJ)/*.o	

