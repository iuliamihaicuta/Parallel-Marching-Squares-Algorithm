build: marching_squares.c
	gcc marching_squares.c helpers.c -o marching_squares -lm -lpthread -Wall -Wextra
clean:
	rm -rf marching_squares