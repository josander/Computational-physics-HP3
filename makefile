
CC = gcc
CFLAGS = -O3
LIBS = -lm

HEADERS = func.h
OBJECTS = H3_main_task2.o func.o
PROGRAM = H3

%.o: %.c $(HEADERS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: $(PROGRAM)

$(PROGRAM): $(OBJECTS)
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)

clean:
	rm -f *.o
	touch *.c

