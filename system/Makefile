CFLAGS = -g

CFLAGS += `gdal-config --cflags gdal`
CFLAGS += `pkg-config --cflags opencv`

LFLAGS = -g
LFLAGS += `gdal-config --libs gdal`
LFLAGS += `pkg-config --libs opencv`
LFLAGS += -L libs -l vl

CPPFLAGS = -std=c++11

BINDIR = bin
OBJDIR = obj

SRCS = $(wildcard *.cpp)
OBJS = $(addprefix $(OBJDIR)/,$(SRCS:.cpp=.o))

DESCRIPTORS_OBJS = $(wildcard descriptors/obj/*.o)

EXE = classify

.PHONY: directories descriptors

build: directories descriptors $(OBJS) copylibs
	g++ $(OBJS) $(DESCRIPTORS_OBJS) $(LFLAGS) -o $(BINDIR)/$(EXE)

$(OBJDIR)/%.o : %.cpp
	g++ -c $(CFLAGS) $(CPPFLAGS) $< -o $@

directories: $(OBJDIR) $(BINDIR)

$(OBJDIR):
	mkdir -p $(OBJDIR)

$(BINDIR):
	mkdir -p $(BINDIR)

descriptors:
	cd descriptors && $(MAKE)

copylibs:
	cp -f libs/libvl.dylib bin/

run:
	./run.sh

clean:
	rm -rf $(OBJDIR)
	rm -rf $(BINDIR)
	cd descriptors && $(MAKE) clean

