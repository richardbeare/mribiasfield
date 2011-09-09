#include "fftw3.h"
#include "itkImage.h"
#include "Debug.h"
#include "ioutils.h"

#include "itkNeighborhoodOperatorImageFilter.h"
#include "itkBackwardDifferenceOperator.h"
#include "itkForwardDifferenceOperator.h"
#include "itkReflectiveBoundaryCondition.h"
#include "itkConstantBoundaryCondition.h"
#include "itkImportImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkDivideByConstantImageFilter.h"

#define NUM_OF_THREADS 1 //set number of threads to number of available processors

#define M_PI_FLOAT 3.14159265358979323846f


#include "itkNeighborhoodAlgorithm.h"

template <class TImage>
void resetBorder(typename TImage::Pointer dest, typename TImage::Pointer source, unsigned dimension, float multiplier)
{
  // copy parts of the boundary from one image to another
  typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<TImage>::FaceListType faceList;
  itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<TImage> fC;
  typename TImage::SizeType radius;
  radius.Fill(0);
  radius[dimension]=1;
  faceList = fC(source, source->GetLargestPossibleRegion(), radius);

  typename itk::NeighborhoodAlgorithm::ImageBoundaryFacesCalculator<TImage>::FaceListType::iterator fit;

  itk::ImageRegionIterator<TImage> o_iter, i_iter;

  fit = faceList.begin();
  ++fit; // first is the body

  for ( ;fit !=faceList.end(); ++fit)
    {
    typename TImage::SizeType sz = fit->GetSize();
    if (sz[dimension]==1)
      {
      // copy this one
      o_iter = itk::ImageRegionIterator<TImage>(dest, *fit);
      i_iter = itk::ImageRegionIterator<TImage>(source, *fit);
      o_iter.GoToBegin();
      i_iter.GoToBegin();

      for (; !o_iter.IsAtEnd(); ++o_iter, ++i_iter)
	{
	o_iter.Set(multiplier * i_iter.Get());
	}

      }
    }

}

#include "itkImageRegionIteratorWithIndex.h"

template <class TImage>
void fourierSolve(typename TImage::Pointer& imgData, typename TImage::Pointer& imgGradX, 
		  typename TImage::Pointer& imgGradY, typename TImage::Pointer& imgGradZ, float dataCost)
{

  unsigned threadCount = itk::MultiThreader::GetGlobalMaximumNumberOfThreads();
  int retVal = fftwf_init_threads();
  INSIST(retVal != 0);
  //CShape imgShape = imgData.Shape();
  //int nodeCount   = imgShape.NodeCount();
  int nodeCount = imgData->GetLargestPossibleRegion().GetNumberOfPixels();
  typename TImage::SizeType sz = imgData->GetLargestPossibleRegion().GetSize();
  int width = sz[0];
  int height= sz[1];
  int nBands = 1;

  INSIST(nodeCount > 0);

  float* fftBuff   = (float*) fftwf_malloc(sizeof(*fftBuff)   * nodeCount);
  INSIST(fftBuff   != NULL);

  typedef typename itk::ImportImageFilter<typename TImage::PixelType, TImage::ImageDimension> ImportType;
  typename ImportType::Pointer importer = ImportType::New();
  importer->SetRegion(imgData->GetLargestPossibleRegion());
  importer->SetImportPointer(fftBuff, (unsigned long)nodeCount, false);

  typename TImage::Pointer fftIm = importer->GetOutput();
  fftIm->Update();
  fftIm->DisconnectPipeline();

  fftIm->CopyInformation(imgData);

  // put the gradient images in a vector
  std::vector<typename TImage::Pointer> GradImVec(TImage::ImageDimension);
  GradImVec[0]=imgGradX;
  GradImVec[1]=imgGradY;
  if (TImage::ImageDimension == 3)
    {
    GradImVec[2]=imgGradZ;
    }

  //compute two 1D lookup tables for computing the DCT of a 2D Laplacian on the fly
  // float* ftLapY = (float*) fftwf_malloc(sizeof(*ftLapY) * height);
  // float* ftLapX = (float*) fftwf_malloc(sizeof(*ftLapX) * width);

  typename TImage::PointType lapoff;
  lapoff.Fill(0);
  lapoff[TImage::ImageDimension - 1] = -2.0 * TImage::ImageDimension;


  std::cout << lapoff << std::endl;
  std::vector<float *> ftLap(TImage::ImageDimension);
  for (unsigned d = 0; d < TImage::ImageDimension; d++)
    {
    ftLap[d]=(float *)fftwf_malloc(sizeof(*(ftLap[d])) * sz[d]);
    INSIST(ftLap[d] != NULL);

    for (unsigned p = 0; p < sz[d]; p++)
      {
      ftLap[d][p] = 2.0f * cos(M_PI_FLOAT * p / (sz[d] - 1)) + lapoff[d];
      }

    }

  // INSIST(ftLapX != NULL);
  // INSIST(ftLapY != NULL);

  // need to extend this to include ftLapZ

  // for(int x = 0; x < width; x++)
  //   {
  //   ftLapX[x] = 2.0f * cos(M_PI_FLOAT * x / (width - 1));
  //   }
  // for(int y = 0; y < height; y++)
  //   {
  //   ftLapY[y] = -4.0f + (2.0f * cos(M_PI_FLOAT * y / (height - 1)));
  //   }

  //Create a DCT-I plan for, which is its own inverse.
  fftwf_plan_with_nthreads(threadCount);
  fftwf_plan fftPlan;	

  //use FFTW_PATIENT when plan can be reused
  if (TImage::ImageDimension == 2)
    {
    fftPlan = fftwf_plan_r2r_2d(height, width, 
				fftBuff, fftBuff, 
				FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE); 
    }
  else if (TImage::ImageDimension == 3)
    {
    // use patient here, because everything is big
    fftPlan = fftwf_plan_r2r_3d(sz[2], sz[1], sz[0],
				fftBuff, fftBuff, 
				FFTW_REDFT00, FFTW_REDFT00, FFTW_REDFT00, FFTW_ESTIMATE); 

    }
  else
    {
    std::cerr << "Unsupported dimension" << std::endl;
    return;
    }

  for(int iChannel = 0; iChannel < nBands; iChannel++)
    {
    printf("Solving channel - %i\n", iChannel); 

    int nodeAddr        = 0;
    int pixelAddr       = iChannel;
    int rightPixelAddr  = nBands + iChannel;
    int topPixelAddr    = (width * nBands) + iChannel;

    float dcSum = 0.0f;

#if 0
    // compute h_hat from u, gx, gy (see equation 28 in the paper), as well as the DC term of u's DCT.
    for(int y = 0; y < height; y++)
      for(int x = 0; x < width;  x++, 
	    nodeAddr++, pixelAddr += nBands, rightPixelAddr += nBands, topPixelAddr += nBands)
	{
	// Compute DC term of u's DCT without computing the whole DCT.
	float dcMult = 1.0f;
	if((x > 0) && (x < width  - 1))
	  dcMult *= 2.0f;
	if((y > 0) && (y < height - 1))
	  dcMult *= 2.0f;
	dcSum += dcMult * imgData->GetBufferPointer()[pixelAddr];

	fftBuff[nodeAddr] = dataCost * imgData->GetBufferPointer()[pixelAddr];

	// Subtract g^x_x and g^y_y, with boundary factor of -2.0 to account for boundary reflections implicit in the DCT
	if((x > 0) && (x < width - 1))
	  fftBuff[nodeAddr] -= (imgGradX->GetBufferPointer()[rightPixelAddr] - imgGradX->GetBufferPointer()[pixelAddr]);
	else
	  fftBuff[nodeAddr] -= (-2.0f * imgGradX->GetBufferPointer()[pixelAddr]);

	if((y > 0) && (y < height - 1))
	  fftBuff[nodeAddr] -= (imgGradY->GetBufferPointer()[topPixelAddr] - imgGradY->GetBufferPointer()[pixelAddr]);
	else
	  fftBuff[nodeAddr] -= (-2.0f * imgGradY->GetBufferPointer()[pixelAddr]);
	}
#else
    // in an ITK form
    typedef itk::ImageRegionIteratorWithIndex<TImage> IndexItType;
    typedef itk::ImageRegionIterator<TImage> ItType;
    IndexItType imBuffIt(imgData, imgData->GetLargestPossibleRegion());
    ItType fftBuffIt(fftIm, fftIm->GetLargestPossibleRegion());

    imBuffIt.GoToBegin();
    fftBuffIt.GoToBegin();


    for (;!imBuffIt.IsAtEnd();++imBuffIt, ++fftBuffIt)
      {
      typename TImage::IndexType here = imBuffIt.GetIndex();
      float dcMult = 1.0f;
      for (unsigned dim = 0; dim < TImage::ImageDimension; dim++)
	{
	if ((here[dim] > 0) && (here[dim] < ((int)sz[dim] - 1)))
	  {
	  dcMult *= 2.0f;
	  }

	}
      typename TImage::PixelType imval = imBuffIt.Get();
      dcSum += dcMult * imval;
      fftBuffIt.Set(dataCost * imval);

      }

    // could compute all of the forward difference images, then set
    // edges to modified values - we are going to reset the
    // boundaries, so use a constant boundary condition.
    std::vector<typename TImage::Pointer> DGradImVec(TImage::ImageDimension);


    typedef typename itk::ForwardDifferenceOperator<typename TImage::PixelType, TImage::ImageDimension> ForwardOper;
    typedef typename itk::NeighborhoodOperatorImageFilter<TImage, TImage> OpFilterType;

    itk::ConstantBoundaryCondition<TImage> cbc;
    typename OpFilterType::Pointer difference = OpFilterType::New();

    difference->OverrideBoundaryCondition(&cbc);

    for (unsigned dim = 0; dim < TImage::ImageDimension; dim++)
      {
      ForwardOper dd;
      dd.SetDirection(dim);
      dd.CreateDirectional();
      
      difference->SetInput(GradImVec[dim]);
      difference->SetOperator(dd);
      DGradImVec[dim]=difference->GetOutput();
      DGradImVec[dim]->Update();
      DGradImVec[dim]->DisconnectPipeline();
      resetBorder<TImage>(DGradImVec[dim], GradImVec[dim], dim, -2.0);
      }

//    writeIm<TImage>(DGradImVec[2], "tmp.mha");

    typedef typename itk::SubtractImageFilter<TImage, TImage, TImage> SubtractType;
    typename SubtractType::Pointer subtract = SubtractType::New();
    for (unsigned dim = 0; dim < TImage::ImageDimension; dim++)
      {
      subtract->InPlaceOn();
      subtract->SetInput1(fftIm);
      subtract->SetInput2(DGradImVec[dim]);
      subtract->Modified();
      subtract->Update();
      fftIm = subtract->GetOutput();
      fftIm->DisconnectPipeline();
      }
 
#endif

    writeIm<TImage>(fftIm, "fftIm.nii.gz");
    //transform h_hat to H_hat by taking the DCT of h_hat
    fftwf_execute(fftPlan);

    
    //compute F_hat using H_hat (see equation 29 in the paper)
    nodeAddr = 0;
    // about the last non dimension independent part, and non ITK part
#if 0
    for(int y = 0; y < height; y++)
      for(int x = 0; x < width;  x++, nodeAddr++)
	{
	//float ftLapResponse = ftLapY[y] + ftLapX[x]; 
	float ftLapResponse = ftLap[0][y] + ftLap[1][x]; 
	fftBuff[nodeAddr] /= (dataCost - ftLapResponse);
	}
#else
    // only dimension independent way is with an  index iterator?
//    typedef itk::ImageRegionIteratorWithIndex<TImage> IndexItType;
    IndexItType ftBuffIt(fftIm, fftIm->GetLargestPossibleRegion());
    
    for (ftBuffIt.GoToBegin(); !ftBuffIt.IsAtEnd(); ++ftBuffIt)
      {
      typename TImage::IndexType here = ftBuffIt.GetIndex();
      float ftLapResponse = 0.0;
      // for (unsigned k=0; k<TImage::ImageDimension; k++)
      // 	{
      // 	ftLapResponse += ftLap[k][here[k]];
      // 	}
      ftLapResponse = ftLap[0][here[0]] + ftLap[1][here[1]] + ftLap[2][here[2]];
      ftBuffIt.Set(ftBuffIt.Get()/(dataCost - ftLapResponse));
      }

#endif
    /**
		 * Set the DC term of the solution to the value computed above (i.e., the DC term of imgData). 
		 * When dataCost = 0 (i.e., there is no data image and the problem becomes pure gradient field integration)
		 * then the DC term  of the solution is undefined. So if you want to control the DC of the solution 
		 * when dataCost = 0 then before calling fourierSolve() set every pixel in 'imgData' to the average value 
		 * you would like the pixels in the solution to have. 
		 */
    fftBuff[0] = dcSum;
    
    //transform F_hat to f_hat by taking the inverse DCT of F_hat
    fftwf_execute(fftPlan);		
    float fftDenom = 4.0f * (width - 1) * (height - 1);
    pixelAddr = iChannel;
#if 0
    for(int iNode = 0; iNode < nodeCount; iNode++, pixelAddr += nBands)
      {
      imgData->GetBufferPointer()[pixelAddr] = fftBuff[iNode] / fftDenom;	
      }
#else
  
    typedef typename itk::DivideByConstantImageFilter<TImage, float, TImage> DivType;
    typename DivType::Pointer divider = DivType::New();
    divider->SetInput(fftIm);
    divider->SetConstant(fftDenom);
    imgData = divider->GetOutput();
    imgData->Update();
    imgData->DisconnectPipeline();

#endif
    }

  fftwf_free(fftBuff);
  // fftwf_free(ftLapX);
  // fftwf_free(ftLapY);
  for (unsigned d=0; d < TImage::ImageDimension; d++)
    {
    fftwf_free(ftLap[d]);
    }
  fftwf_destroy_plan(fftPlan);
  fftwf_cleanup_threads();

  printf("\tDone.\n");
}


////////////////////////////////////////////////////////
int main(int argc, char * argv[])
{
  typedef float PixType;
  const unsigned dim = 3;

  typedef itk::Image<PixType, dim> ImType;

  ImType::Pointer im, g1, g2, g3;
  itk::MultiThreader::SetGlobalMaximumNumberOfThreads(1);

  im = readIm<ImType>(argv[1]);

  typedef itk::NeighborhoodOperatorImageFilter<ImType, ImType> OpFilterType;
  OpFilterType::Pointer difference = OpFilterType::New();
  
  itk::ReflectiveBoundaryCondition<ImType> rbc;
  rbc.SetReflectionType(3);

  typedef itk::BackwardDifferenceOperator<PixType, dim> BackOper;
  typedef itk::ForwardDifferenceOperator<PixType, dim> ForwardOper;
  BackOper bdX, bdY, bdZ;
  ForwardOper fdY;

  bdX.SetDirection(0);
  bdX.CreateDirectional();
  difference->SetOperator(bdX);
  difference->OverrideBoundaryCondition(&rbc);

  difference->SetInput(im);

  g1 = difference->GetOutput();
  g1->Update();
  g1->DisconnectPipeline();


  bdY.SetDirection(1);
  bdY.CreateDirectional();
  difference->SetOperator(bdY);
  difference->OverrideBoundaryCondition(&rbc);

  g2 = difference->GetOutput();
  g2->Update();
  g2->DisconnectPipeline();


  bdZ.SetDirection(2);
  bdZ.CreateDirectional();
  // // std::cout << bdZ << std::endl;
  // // std::cout << bdY << std::endl;

  difference->SetOperator(bdZ);
  difference->OverrideBoundaryCondition(&rbc);

  g3 = difference->GetOutput();
  g3->Update();
  g3->DisconnectPipeline();
  
  // writeIm<ImType>(g1, "dX.nii.gz");
  // writeIm<ImType>(g2, "dY.nii.gz");
  // writeIm<ImType>(g3, "dZ.nii.gz");


  im->FillBuffer(100);
  fourierSolve<ImType>(im, g1, g2, g3, 0);
  writeIm<ImType>(im, "output.mha");

  return(EXIT_SUCCESS);
}

