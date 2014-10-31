CFLAGS = 

CFLAGS += `pkg-config --cflags gdal`
CFLAGS += `pkg-config --cflags opencv`

LFLAGS =
LFLAGS += `pkg-config --libs gdal`
LFLAGS += `pkg-config --libs opencv`

BINDIR = bin
OBJDIR = obj

SRCS = $(wildcard *.cpp)
OBJS = $(addprefix $(OBJDIR)/,$(SRCS:.cpp=.o))

EXE = classify

.PHONY: directories

build: directories $(OBJS)
	g++ $(OBJS) $(LFLAGS) -o $(BINDIR)/$(EXE)

$(OBJDIR)/%.o : %.cpp
	g++ -c $(CFLAGS) $< -o $@

directories: $(OBJDIR) $(BINDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

run:
	./$(BINDIR)/$(EXE)

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
