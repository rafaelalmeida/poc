CFLAGS = 

CFLAGS += `pkg-config --cflags gdal`
CFLAGS += `pkg-config --cflags opencv`

LFLAGS =
LFLAGS += `pkg-config --libs gdal`
LFLAGS += `pkg-config --libs opencv`

EXE = classify

BINDIR = bin
OBJDIR = obj

SRCS=$(wildcard *.cpp)
OBJS=$(SRCS:.cpp=.o)

build: $(OBJS)
	mkdir -p $(BINDIR)
	g++ $(OBJDIR)/*.o $(LFLAGS) -o $(BINDIR)/$(EXE)

%.o: $(SRCS)
	mkdir -p $(OBJDIR)
	g++ -c $(CFLAGS) -o $(OBJDIR)/$@ $^

run:
	./$(BINDIR)/$(EXE)

clean:
	rm -rf $(BINDIR)
	rm -rf $(OBJDIR)
	rm -f $(OBJS)