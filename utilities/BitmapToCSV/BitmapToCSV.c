#include "BitmapToCSV.h"


int main(int argc, char *argv[])
{
	int opt;
	char *inputFileName;
	char *outputFileName;
	int inputs = 0;
	
	FILE *inputFile;
	FILE *outputFile;
	
	Pixel *pixelArray;
	int totalPixels;
	
	//Get command line arguments
	while((opt = getopt(argc, argv, "i:o:")) != -1)
	{
		switch(opt)
		{
			case 'i':
				inputFileName = optarg;
				inputs++;
			break;
			
			case 'o':
				outputFileName = optarg;
				inputs++;
			break;
		}
	}
	//Check that input and output files were specified
	if(inputs < 2)
	{
		printf("Usage: ./BitmapToCSV -i INPUT_FILE -o OUTPUT_FILE\n");
		return 1;
	}
	
	//Open input file
	inputFile = fopen(inputFileName, "r");
	//Get pixel array
	getPixelData(inputFile, pixelArray, &totalPixels);
	//Open output file
	
	//Convert pixel data to csv and output to output file
	
}

void getPixelData(FILE *bitmapFile, Pixel *pixelArray, int *pixelsRead)
{
	uint32_t width;
	uint32_t height;
	
	//Read entire file into buffer
	fseek(bitmapFile, 0, SEEK_END);
	long length = ftell(bitmapFile);
	char *buffer = malloc(length);	
	rewind(bitmapFile);
	fread(buffer, length, 1, bitmapFile);
	
	//Get width (at position 0x12)
	char *tempPointer = buffer + 0x12;
	//Cast down to 32 bit signed int
	width = (uint32_t)(*tempPointer);
	printf("width: %d\n", width);
	
	//Get height (at position 0x16)
	tempPointer = buffer + 0x16;
	height = (uint32_t)(*tempPointer);
	printf("height: %d\n", height);
	
	//Allocate memory for pixel array
	*pixelsRead = width * height;
	pixelArray = (Pixel *)malloc(width*height*sizeof(Pixel));
	
	char *currentByte = buffer + 0x36;
	int offset = 0;
	for(int y = 0, i = 0; y < height; y++)
	{
		for(int x = 0; x < width; x++, i++)
		{
			pixelArray[i].x = x;
			pixelArray[i].y = y;
			//Set color values
			pixelArray[i].r = *(currentByte + offset);
			pixelArray[i].g = *(currentByte + 1 + offset);
			pixelArray[i].b = *(currentByte + 2 + offset);
			
			printf("x:%d, y:%d, r:%dg:%d,b:%d\n", x, y, (uint8_t)pixelArray[i].r, (uint8_t)pixelArray[i].g, (uint8_t)pixelArray[i].b);
			offset += 3;
		}
		//Account for alignment / padding
		//Each row must be a multiple of 4
		//Probably a better way to do this but it works
		while(offset % 4 != 0)
			offset++;
	}
	
	free(buffer);
}
