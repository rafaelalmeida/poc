CFLAGS = -O2
LFLAGS =

BINDIR = bin
OBJDIR = obj

SRCS = $(wildcard *.c)
OBJS = $(addprefix $(OBJDIR)/,$(SRCS:.c=.o))

.PHONY: directories

objs: directories $(OBJS)

$(OBJDIR)/%.o : %.c
	clang -c $(CFLAGS) $< -o $@

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
