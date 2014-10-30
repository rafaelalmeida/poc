CFLAGS = 

CFLAGS += `pkg-config --cflags gdal`
CFLAGS += `pkg-config --cflags opencv`

LFLAGS =
LFLAGS += `pkg-config --libs gdal`
LFLAGS += `pkg-config --libs opencv`

EXE = classify

BINDIR = bin
OBJDIR = obj

build: gdal_driver main
	mkdir -p $(BINDIR)
	g++ $(OBJDIR)/*.o $(LFLAGS) -o $(BINDIR)/$(EXE)

main: main.cpp
	mkdir -p $(OBJDIR)
	g++ -o $(OBJDIR)/main.o -c $(CFLAGS) main.cpp

gdal_driver: gdal_driver.cpp
	mkdir -p $(OBJDIR)
	g++ -o $(OBJDIR)/gdal_driver.o -c $(CFLAGS) gdal_driver.cpp

run:
	./$(BINDIR)/$(EXE)

clean:
	rm -rf $(BINDIR)
	rm -rf $(OBJDIR)
	rm -f $(OBJS)