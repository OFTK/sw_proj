CC = gcc
CC_CMD = $(CC) -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c -lm -o spkmeans

makec:
	$(CC_CMD) 


debug: CC_CMD += -g
debug: makec

all: makec

clean:
	rm spkmeans
