#include <iostream>
#include <armadillo>
#include <math.h>   // This gives INFINITY
#include <typeinfo>
#include <stack>
#include <ctime>



#ifdef __cplusplus 
extern "C" {
#endif
    /* Declarations of this file */

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <signal.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "adcdriver_host.h"
#include "spidriver_host.h"
#ifdef __cplusplus
}
#endif


using namespace std;
using namespace arma;

const float pi = 3.1415926535897932384;
#define NUM_FRAMES 100 //100
#define NUM_PTS 1000 //1024
#define MAX_RANK 20
#define ZEROTHRESH 1.0e-7
fvec v_sample[NUM_PTS];
float u[NUM_PTS];
char dummy[8];

int main(int argc, char** argv){
	
	// Buffers for tx and rx data from A/D registers
	uint32_t tx_buf[3];
	uint32_t rx_buf[4];
	// Initialize A/D converter
	adc_config();
	// Check the A/D is alive by reading from its config reg.
	printf("--------------------------------------------------\n");
	printf("About to read A/D config register\n");
	printf("Hit return when ready -->\n");
	fgets (dummy, 8, stdin);
	rx_buf[0] = adc_get_id_reg();
	printf("Read ID reg.  Received ID = 0x%08x\n", rx_buf[0]);
	
	// Set sample rate
	printf("--------------------------------------------------\n");
	printf("Set sample rate and set channel 0\n");
	adc_set_samplerate(SAMP_RATE_31250);
	adc_set_chan0();
	
	printf("--------------------------------------------------\n");
	printf("Now read training set: %d data frames\n", NUM_FRAMES);
	printf("Hit return when ready -->\n");
	fgets (dummy, 8, stdin);
	
	// Do multiple read
	printf("--------------------------------------------------\n");
	printf("Now try to do multiple read\n");
	printf("Hit return when ready -->\n");
	fgets (dummy, 8, stdin);
	
	/*
	// This is manuyally created signal
	int Nx = NUM_PTS + NUM_FRAMES - 1;
	
	int Tmin = 0; int Tmax = 3;
	float rbound = Tmax*(Nx-1)/Nx;
	fvec t = linspace<fvec>(Tmin,rbound,Nx);
	
	float dt;
	dt = t(1) - t(0);
	float fsamp = 1/dt;
	
	float f0;
	f0 = 2.0;
	fvec v_sample = sin(2*pi*f0*t);
	*/
	printf("--------------------------------------------------\n");
	printf("Preparations are done and now we are going to the main part of the algorithm\n");
	printf("Hit return when ready -->\n");
	fgets (dummy, 8, stdin);
	
	//adc_read_multiple(NUM_PTS, u); 
	
	
	fmat Yf;
	fvec ypsd;
	// Loop over the signals data
	for(int i = 0; i < NUM_FRAMES; i++){
		adc_read_multiple(NUM_PTS, u); 
		
		fvec v(NUM_PTS);
		//v = v_sample(span(i:NUM_PTS+i-1));
		for(int j = 0; j < NUM_PTS; j++){
			v(j) = u[j];
		}
		
		// Compute mean
		float v_mean = mean(v);
		// Substruct mean
		fvec y = v - v_mean;
		
		// Do the fft
		cx_fvec yspec;
		yspec = fft(y);
		// Compute the power spectral density
		//fvec ypsd;
		ypsd = abs(yspec);
		ypsd = square(ypsd);
		
		// Normalize and threshold spectrum
		float tmp1 = sum(ypsd);
		float pow_learn_thresh = tmp1/1000.f;
		
		ypsd = ypsd/tmp1;
		//fvec idx = find(ypsd < ZEROTHRESH);
		//for(int j = 0; j<size(idx))
		for(int j = 0; j < size(ypsd,0); j++){
			if(ypsd(j) < ZEROTHRESH){
				ypsd(j) = 0.0f;
			}
		}
		// Create the matrix Yf
		Yf = join_rows(Yf,ypsd);
	}
	// Do svd
	fmat U, V; fvec s;
	svd(U,s,V,Yf);
	// Get effective rank of M0, then pull out vectors
	float pr = sum(s);
	fvec s_over_sum = s/pr;
	float H = 0;
	for(int i = 0; i < NUM_FRAMES; i++){
		if(s_over_sum(i) > 0){
			float entry = s_over_sum(i);
			H = H + entry * log(entry);
		}
	}
	int r = floor(exp(-H));
	cout << "r = " << r << "\n\n";
	//r = 1;
	// If the dimensionality beyond the maximum, set dimen = r
	if(r > MAX_RANK){
		r = MAX_RANK;
	}
	
	// Extract
	fmat Ur = U.cols(0,r-1); // The index starts at 0, so minus down by 1
	// Create projection
	// Nr = I - Ur * ( (Ur' * Ur) \ Ur' );
	fmat Urtrans = Ur.t(); // Transpose, checked
	fmat I = eye<fmat>(size(Ur,0),size(Ur,0));
	fmat A = inv(Urtrans*Ur);
	fmat Nr = I - Ur * (A*Urtrans);
	
	//cout << "ypsd = " << ypsd << endl;
	//cout << "norm(ypsd) = " << norm(ypsd,2) << endl;
	
	
	fvec Nrmypsd = Nr*ypsd;
	fvec newypsd = abs(ypsd);
	float n1 = 0; float n2 = 0;
	for(int j = 0; j<size(ypsd,0); j++){
		n1 += Nrmypsd(j) * Nrmypsd(j);
		n2 += newypsd(j) * newypsd(j);
	}
	n1 = sqrt(n1);
	n2 = sqrt(n2);
	
	
	
	
	
	// Compute the norm
	float yn = norm(Nr*ypsd)/norm(ypsd);
	yn = n1/n2;
	// Use K = 2.8 for example
    	float K = 2.8;
    	// Threshold
    	float Tr = K*yn;
    
    	cout << "The threshold Tr is " << Tr << endl;
    
    
    	// --------------------------------------
    	cout << " Now do the operating phase\n";
    	// Take a segment of the sine wave in the training phase
    	// Suppose v is the wave segment
    	//fvec t1 = linspace<fvec>(Tmin,rbound,NUM_PTS);
    	/*
    	fvec t1(NUM_PTS);
    	for(int i = 0; i<NUM_PTS; i++){
    		t1(i) = t(i);
    	}
    	int f1 = 2;
    	fvec v = sin(2*pi*f1*t1);
    	*/
    
    	adc_read_multiple(NUM_PTS, u); 
    	fvec v(NUM_PTS);
	for(int i = 0; i < NUM_PTS; i++){
		v(i) = u[i];
	}
    
    
    	fvec x = v - mean(v);
    	// Compute fft
    	cx_fvec xspec = fft(x);
    	fvec xpsd = abs(xspec);
    	xpsd = square(xpsd);
    	float tmp1 = sum(xpsd);
    	xpsd = xpsd/tmp1;
    	// Normalize and threshold spectrum
    	for(int j = 0; j< size(xpsd,0); j++){
		if(v(j) < ZEROTHRESH){
			v(j) = 0.0f;
		}
	}
	
	// Project onto nullspace
	fvec xproj = Nr * xpsd;
	// Get the norm of xproj
	float xprojNorm = norm(xproj,2);
	// Get the norm of xproj
	float xpsdNorm = norm(xpsd,2);
	float ratioNorm = xprojNorm/xpsdNorm;
	cout << "The ratio of two norms is " << ratioNorm;
	cout << " and Tr = " << Tr << endl;
	if(ratioNorm > Tr){
		cout << "Anomaly detected\n";
	}
	else{
		cout << "Seems like it is normal\n";
	}

	
}
