#ifndef __sprMRnormalization_hxx
#define __sprMRnormalization_hxx

#include "sprMRnormalization.h"
#include "itkObjectFactory.h"

#include "itkImageToHistogramFilter.h"
#include "itkShiftScaleImageFilter.h"
#include "itkThresholdImageFilter.h"

#include "itkImageRegionIterator.h"

namespace spr
{

template< class TImage>
void mrNormalization< TImage>
::GenerateData()
{
  typename TImage::ConstPointer input = this->GetInput();
  typename TImage::Pointer output = this->GetOutput();
  this->AllocateOutputs();

  typedef itk::Statistics::ImageToHistogramFilter< TImage > HistogramFilterType;
  typename HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();
  typename HistogramFilterType::HistogramSizeType numBins( 1 );
  numBins[0] = 100;
  histogramFilter->SetHistogramSize( numBins );
  histogramFilter->SetMarginalScale( 10.0 );
  histogramFilter->SetAutoMinimumMaximum(true);
  histogramFilter->SetInput( input );
  histogramFilter->Update();

  typename HistogramFilterType::HistogramType * histogram = histogramFilter->GetOutput();
  typename TImage::PixelType min,max;
  min = histogram->Quantile(0,tl);
  max = histogram->Quantile(0,tu);

  itk::ImageRegionConstIterator<TImage> itInput(input, input->GetLargestPossibleRegion());
  itk::ImageRegionIterator<TImage> itOutput(output, output->GetLargestPossibleRegion());

  typename TImage::PixelType tmp;
  for (itInput.GoToBegin(), itOutput.GoToBegin(); !itInput.IsAtEnd(); ++itInput, ++itOutput)
  {
	  tmp = (itInput.Get() - min) / (max - min);
	  tmp = (tmp > 1) ? 1 : tmp;
	  tmp = (tmp < 0) ? 0 : tmp;
	  itOutput.Set(tmp);
  }
}

}// end namespace


#endif
