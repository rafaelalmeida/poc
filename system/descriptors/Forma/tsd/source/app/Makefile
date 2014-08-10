INC=../include
LIB=../lib
BIN=../bin
FLAGS= -O2 -Wall

$@.c: $(LIB)/libift.a $@.c 
	cd ../; make; cd app
	gcc $(FLAGS) $@.c -o $(BIN)/$@ -I$(INC) -L$(LIB) -lift -lm

clean:
	rm $(BIN)/*

