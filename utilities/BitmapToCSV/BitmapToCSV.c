#include "BitmapToCSV.h"

/*
Utility for creating test datasets from the pixel data of bitmap files.
Any pixel that is not white will be added to the data set using its position.
Bitmap file must be 24 bit color with color space information not included
*/

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
	while((opt = getopt(argc, argv, "i:o:h")) != -1)
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
			
			case 'h':
				printf("Must include input file and output file\n \t-i: INPUT_FILE_NAME\n\t-o: OUTPUT_FILE_NAME\n");
				return 1;
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
	getPixelData(inputFile, &pixelArray, &totalPixels);
	//Open output file
	outputFile = fopen(outputFileName, "w");
	//Convert pixel data to csv and output to output file
	outputCSV(outputFile, pixelArray, totalPixels);
	
	fclose(inputFile);
	fclose(outputFile);
}

void outputCSV(FILE *outputFile, Pixel *pixelArray, int pixelArraySize)
{
	for(int i = 0; i < pixelArraySize; i++)
	{
		//Check if the current pixel is not white
		if(((uint8_t)pixelArray[i].r) != 0xFF || ((uint8_t)pixelArray[i].g) != 0xFF || ((uint8_t)pixelArray[i].b) != 0xFF)
		{
			//Output to csv file as point
			fprintf(outputFile, "%d,%d\n", pixelArray[i].x, pixelArray[i].y);
			//printf("outputing\n");
		}
	}
}

void getPixelData(FILE *bitmapFile, Pixel **pixelArray, int *pixelsRead)
{
	int32_t width;
	int32_t height;
	
	//Read entire file into buffer
	fseek(bitmapFile, 0, SEEK_END);
	long length = ftell(bitmapFile);
	char *buffer = malloc(length);	
	rewind(bitmapFile);
	fread(buffer, length, 1, bitmapFile);
	
	//Get width (at position 0x12)
	char *tempPointer = buffer + 0x12;
	//Cast down to 32 bit signed int
	width = *((int32_t*)tempPointer);
	printf("width: %d\n", width);
	
	//Get height (at position 0x16)
	tempPointer = buffer + 0x16;
	height = *((int32_t*)tempPointer);
	printf("height: %d\n", height);
	
	//Allocate memory for pixel array
	*pixelsRead = width * height;
	*pixelArray = (Pixel *)malloc(width*height*sizeof(Pixel));
	
	char *currentByte = buffer + 0x36;
	int offset = 0;
	for(int y = 0, i = 0; y < height; y++)
	{
		for(int x = 0; x < width; x++, i++)
		{
			(*pixelArray)[i].x = x;
			(*pixelArray)[i].y = y;
			//Set color values
			(*pixelArray)[i].r = *(currentByte + offset);
			(*pixelArray)[i].g = *(currentByte + 1 + offset);
			(*pixelArray)[i].b = *(currentByte + 2 + offset);
			
			printf("x:%d, y:%d, r:%dg:%d,b:%d\n", x, y, (uint8_t)(*pixelArray)[i].r, (uint8_t)(*pixelArray)[i].g, (uint8_t)(*pixelArray)[i].b);
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
