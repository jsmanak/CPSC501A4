#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>

using namespace std;

char ChunkID[4];
int ChunkSize;
char format[2];
char subChunk1ID[4];
int subChunk1Size;
int16_t audioFormat;
int16_t numChannels;
int sampleRate;	
int byteRate;
int16_t blockAlign;
int16_t bitsPerSample;
char subChunk2ID[4];
int subChunk2Size;
short* fileData;

int size;

float *readWAV(char *filename, float *signal);
void convolve(float x[], int N, float h[], int M, float y[], int P);
void writeWAV(char *filename, float *signal, int signalSize);


int main(int argc, char **argv)
{

	if (argc < 4) {
		printf("Usage: convolve.exe <input wav> <impulse wav> <output wav>\n");
		return -1;
	}
	
	char *inputName = argv[1];
	char *impulseName = argv[2];
	char *outputName = argv[3];
	
	float *inputSignal = NULL;
	float *IRSignal = NULL;
	float *outputSignal = NULL;
	int inputSignalSize;
	int IRSignalSize;
	int outputSignalSize;
	
	inputSignal = readWAV(inputName, inputSignal);
	inputSignalSize = size;
	printf("\n");
	
	IRSignal = readWAV(impulseName, IRSignal);
	IRSignalSize = size;
	
	outputSignalSize = inputSignalSize + IRSignalSize - 1;
	outputSignal = new float[outputSignalSize];
	
	printf("\n");
	printf("Convolving...\n");
	
	convolve(inputSignal, inputSignalSize, IRSignal, IRSignalSize, outputSignal, outputSignalSize);
	printf("Finished convolution!\n");

	writeWAV(outputName, outputSignal, outputSignalSize);
	
	printf("Complete\n");
	
	return 0;
	
}



	
float *readWAV (char *filename, float *signal) 
{
	
	ifstream inputfile(filename, ios::in | ios::binary);
	
	inputfile.seekg(ios::beg);
	
	inputfile.read(ChunkID, 4);
	inputfile.read((char*) &ChunkSize, 4);
	inputfile.read(format, 4);
	inputfile.read(subChunk1ID, 4);
	inputfile.read((char*) &subChunk1Size, 4);
	inputfile.read((char*) &audioFormat, 2);
	inputfile.read((char*) &numChannels, 2);
	inputfile.read((char*) &sampleRate, 4);
	inputfile.read((char*) &byteRate, 4);
	inputfile.read((char*) &blockAlign, 2);
	inputfile.read((char*) &bitsPerSample, 2);
	if (subChunk1Size == 18) {
		inputfile.seekg(2, ios::cur);
	}
	inputfile.read(subChunk2ID, 4);
	inputfile.read((char*)&subChunk2Size, 4);
	size = subChunk2Size / 2;
	
	short *data = new short[size];
	for (int i = 0; i < size; i++) {
		inputfile.read((char *) &data[i], 2);
	}
	
	inputfile.close();
	printf("Input file %s read!\n", filename);
	
	short sample;
	signal = new float[size];
	printf("Size: %d\n", size);
	for (int i = 0; i < size; i++) {
		sample = data[i];
		signal[i] = (sample * 1.0) / (pow(2.0, 15.0) - 1);
		if (signal[i] < -1.0)
			signal[i] = -1.0;
	}
	
	printf("Input file converted to 1 to -1 range!\n");
	return signal;	
}


//Code from tutorial- had to change double back to float
void convolve(float x[], int N, float h[], int M, float y[], int P)
{
	int n, m;
	
	if (P != (N + M - 1)) {
		printf("Output signal vector is the wrong size\n");
		printf("It is %-d, but should be %-d\n", P, (N + M - 1));
		printf("Aborting convolution\n");
		return;
	}
	
	for (n = 0; n < P; n++)
		y[n] = 0.0;
		
	for (n = 0; n < N; n++) {
		for (m = 0; m < M; m++)
			y[n+m] += x[n] * h[m];
	}
	
	printf("Signal convolved!\n");
}



void writeWAV(char *filename, float *signal, int signalSize)
{

	ofstream outputfile(filename, ios::out | ios::binary);
	
	// File corrupted without hardcoded values
	char *ChunkID = "RIFF";
	char *format = "WAVE";
	// PCM = 18 was unnecessary
	subChunk1Size = 16;
	
	subChunk2Size = numChannels * signalSize * (bitsPerSample / 8);
	ChunkSize = subChunk2Size + 36;
	
	outputfile.write(ChunkID, 4);
	outputfile.write((char*) &ChunkSize, 4);
	outputfile.write(format, 4);
	outputfile.write(subChunk1ID, 4);
	outputfile.write((char*) &subChunk1Size, 4);
	outputfile.write((char*) &audioFormat, 2);
	outputfile.write((char*) &numChannels, 2);
	outputfile.write((char*) &sampleRate, 4);
	outputfile.write((char*) &byteRate, 4);
	outputfile.write((char*) &blockAlign, 2);
	outputfile.write((char*) &bitsPerSample, 2);
	outputfile.write(subChunk2ID, 4);
	outputfile.write((char*)&subChunk2Size, 4);
	
	int16_t sample;
	
	// converting float to int between -2^15 to 2^15 - 1
	for(int i = 0; i < signalSize; i++)
	{
		sample = (int16_t)(signal[i] * (pow(2.0, 15.0) - 1));
		outputfile.write((char*)&sample, 2);
	}
	outputfile.close();
}

// This code was created with help from resources by Zahra Sahaf, Leonard Manzara, and Charles Cote
