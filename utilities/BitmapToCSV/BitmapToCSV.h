#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>


typedef struct
{
	int x;
	int y;
	char r;
	char g; 
	char b;
} Pixel;

void getPixelData(FILE *bitmapFile, Pixel **pixelArray, int *pixelsRead);
void outputCSV(FILE *outputFile, Pixel *pixelArray, int pixelArraySize);
