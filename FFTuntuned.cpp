#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
#include <math.h>
#include <fstream>
#include <iostream>

#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

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
void fft(float data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    float wtemp, wr, wpr, wpi, wi, theta;
    float tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
		if (j > i) {
			SWAP(data[j], data[i]);
			SWAP(data[j+1], data[i+1]);
		}
		m = nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
    }

    mmax = 2;
    while (n > mmax) {
		istep = mmax << 1;
		theta = isign * (6.28318530717959 / mmax);
		wtemp = sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = sin(theta);
		wr = 1.0;
		wi = 0.0;
		for (m = 1; m < mmax; m += 2) {
			for (i = m; i <= n; i += istep) {
				j = i + mmax;
				tempr = wr * data[j] - wi * data[j+1];
				tempi = wr * data[j+1] + wi * data[j];
				data[j] = data[i] - tempr;
				data[j+1] = data[i+1] - tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr = (wtemp = wr) * wpr - wi * wpi + wr;
			wi = wi * wpr + wtemp * wpi + wi;
		}
		mmax = istep;
    }
}

// From notes 12-03-2013
void FFTScaling (float signal[], int N)
{
	int k;
	int i;
	for (k = 0, i = 0; k < N; k++, i+=2) {
		signal[i] /= (float)N;
		signal[i+1] /= (float)N;
	}
}


// Uses overlap-add method of FFT
void convolve(float x[], int N, float h[], int M, float y[], int P) 
{
	int paddedArraySize = 1;
	int i = 0;
	
	// For FFT, need array size of a power of 2
	while (paddedArraySize < P) {
		paddedArraySize *= 2;
	}
	
	// enlarging x, h, and Y and padding with zeros (imaginary)
	float *paddedInput = new float[2 * paddedArraySize];
	for (i = 0; i < (N * 2); i+=2) {
		paddedInput[i] = x[i/2];
		paddedInput[i+1] = 0;
	}
	for (; i < paddedArraySize; i++) {
		paddedInput[i] = 0;
	}
	
	float *paddedIR = new float[2 * paddedArraySize];
	for (i = 0; i < (M * 2); i+=2) {
		paddedIR[i] = h[i/2];
		paddedIR[i+1] = 0;
	}
	for (; i < paddedArraySize; i++) {
		paddedIR[i] = 0;
	}
	
	float *paddedOutput = new float[2 * paddedArraySize];
	for (i = 0; i < paddedArraySize; i++) {
		paddedOutput[i] = 0;
	}
		
	// Convert x and h and X and H
	fft((paddedInput - 1), paddedArraySize, 1);
	fft((paddedIR - 1), paddedArraySize, 1);

	// Complex multiplication to get Y values
	for (i = 0; i < (paddedArraySize * 2); i+=2) {
		paddedOutput[i] = (paddedInput[i] * paddedIR[i]) - (paddedInput[i+1] * paddedIR[i+1]);
		paddedOutput[i+1] = (paddedInput[i+1] * paddedIR[i]) + (paddedInput[i] * paddedIR[i+1]);
	}

	// IFFT on Y -> y
	fft((paddedOutput - 1), paddedArraySize, -1);
	
	// FFT scaling
	FFTScaling(paddedOutput, paddedArraySize);
	
	// removing padding
	for (i = 0; i < P; i++) {
		y[i] = paddedOutput[i*2];
	}
}



void writeWAV(char *filename, float *signal, int signalSize)
{

	ofstream outputfile(filename, ios::out | ios::binary);
	
	char *ChunkID = "RIFF";
	char *format = "WAVE";
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
