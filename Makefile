# define objects files
OBJECTS = \
shooting_auto.o

# define compiler and its flags
CC = g++
CFLAGS = -Wall -O -g

shooting: $(OBJECTS)
				$(CC) $(CFLAGS) -o shooting $(OBJECTS)

%.o : %.c
				$(CC) $(CFLAGS) -o $ <

# remove all objects files and program executables
clean:
				rm -f *.o *.O shooting
