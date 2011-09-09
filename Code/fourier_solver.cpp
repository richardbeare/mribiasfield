#include "fftw3.h"
#include "Image.h"
#define NUM_OF_THREADS 1 //set number of threads to number of available processors

/*
 * The following 100 line function is our entire FFT based solver for the screened Poisson equation.
 * This function won't compile without our Image class (not included but you could replace it with your own) and the fftw library (http://www.fftw.org/). 
 */


void fourierSolve(CFloatImage& imgData, CFloatImage& imgGradX, CFloatImage& imgGradY, float dataCost)
{
	int retVal = fftwf_init_threads();
	INSIST(retVal != 0);

	CShape imgShape = imgData.Shape();
	int nodeCount   = imgShape.NodeCount();
	INSIST(nodeCount > 0);

	float* fftBuff   = (float*) fftwf_malloc(sizeof(*fftBuff)   * nodeCount);
	INSIST(fftBuff   != NULL);

	//compute two 1D lookup tables for computing the DCT of a 2D Laplacian on the fly
	float* ftLapY = (float*) fftwf_malloc(sizeof(*ftLapY) * imgShape.height);
	float* ftLapX = (float*) fftwf_malloc(sizeof(*ftLapX) * imgShape.width);
	INSIST(ftLapX != NULL);
	INSIST(ftLapY != NULL);

	for(int x = 0; x < imgShape.width; x++)
	{
		ftLapX[x] = 2.0f * cos(M_PI_FLOAT * x / (imgShape.width - 1));
	}
	for(int y = 0; y < imgShape.height; y++)
	{
		ftLapY[y] = -4.0f + (2.0f * cos(M_PI_FLOAT * y / (imgShape.height - 1)));
	}

	//Create a DCT-I plan for, which is its own inverse.
	fftwf_plan_with_nthreads(NUM_OF_THREADS);
	fftwf_plan fftPlan;	
	fftPlan = fftwf_plan_r2r_2d(imgShape.height, imgShape.width, 
				    fftBuff, fftBuff, 
				    FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE); //use FFTW_PATIENT when plan can be reused

	for(int iChannel = 0; iChannel < imgShape.nBands; iChannel++)
	{
		printf("Solving channel - %i\n", iChannel); 

		int nodeAddr        = 0;
		int pixelAddr       = iChannel;
		int rightPixelAddr  = imgShape.nBands + iChannel;
		int topPixelAddr    = (imgShape.width * imgShape.nBands) + iChannel;

		float dcSum = 0.0f;

		// compute h_hat from u, gx, gy (see equation 28 in the paper), as well as the DC term of u's DCT.
		for(int y = 0; y < imgShape.height; y++)
		for(int x = 0; x < imgShape.width;  x++, 
			nodeAddr++, pixelAddr += imgShape.nBands, rightPixelAddr += imgShape.nBands, topPixelAddr += imgShape.nBands)
		{
			// Compute DC term of u's DCT without computing the whole DCT.
			float dcMult = 1.0f;
			if((x > 0) && (x < imgShape.width  - 1))
				dcMult *= 2.0f;
			if((y > 0) && (y < imgShape.height - 1))
				dcMult *= 2.0f;
			dcSum += dcMult * imgData[pixelAddr];

			fftBuff[nodeAddr] = dataCost * imgData[pixelAddr];			

			// Subtract g^x_x and g^y_y, with boundary factor of -2.0 to account for boundary reflections implicit in the DCT
			if((x > 0) && (x < imgShape.width - 1))
				fftBuff[nodeAddr] -= (imgGradX[rightPixelAddr] - imgGradX[pixelAddr]);
			else
				fftBuff[nodeAddr] -= (-2.0f * imgGradX[pixelAddr]);

			if((y > 0) && (y < imgShape.height - 1))
				fftBuff[nodeAddr] -= (imgGradY[topPixelAddr] - imgGradY[pixelAddr]);
			else
				fftBuff[nodeAddr] -= (-2.0f * imgGradY[pixelAddr]);
		}

		//transform h_hat to H_hat by taking the DCT of h_hat
		fftwf_execute(fftPlan);

		//compute F_hat using H_hat (see equation 29 in the paper)
		nodeAddr = 0;
		for(int y = 0; y < imgShape.height; y++)
		for(int x = 0; x < imgShape.width;  x++, nodeAddr++)
		{
			float ftLapResponse = ftLapY[y] + ftLapX[x]; 
			fftBuff[nodeAddr] /= (dataCost - ftLapResponse);
		}

		/*
		 * Set the DC term of the solution to the value computed above (i.e., the DC term of imgData). 
		 * When dataCost = 0 (i.e., there is no data image and the problem becomes pure gradient field integration)
		 * then the DC term  of the solution is undefined. So if you want to control the DC of the solution 
		 * when dataCost = 0 then before calling fourierSolve() set every pixel in 'imgData' to the average value 
		 * you would like the pixels in the solution to have. 
		 */
		fftBuff[0] = dcSum;

		//transform F_hat to f_hat by taking the inverse DCT of F_hat
		fftwf_execute(fftPlan);		
		float fftDenom = 4.0f * (imgShape.width - 1) * (imgShape.height - 1);
		pixelAddr = iChannel;
		for(int iNode = 0; iNode < nodeCount; iNode++, pixelAddr += imgShape.nBands)
		{
			imgData[pixelAddr] = fftBuff[iNode] / fftDenom;	
		}
	}

	fftwf_free(fftBuff);
	fftwf_free(ftLapX);
	fftwf_free(ftLapY);
	fftwf_destroy_plan(fftPlan);
	fftwf_cleanup_threads();

	printf("\tDone.\n");
}
