CFLAGS = 

CFLAGS += `pkg-config --libs --cflags gdal`
CFLAGS += `pkg-config --libs --cflags opencv`

build:
	g++ main.cpp $(CFLAGS) -o pipeline

run:
	./pipeline

clean:
	rm -f pipeline