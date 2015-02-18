#ifndef __sprBlockImageFilter_hxx
#define __sprBlockImageFilter_hxx

#include "sprBlockImageFilter.h"
#include "itkObjectFactory.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageToListSampleAdaptor.h"
#include "itkImage.h"
#include "itkMinimumMaximumImageCalculator.h"
#include "itkImageRegionIteratorWithIndex.h"


namespace spr
{

template< class TImage>
void BlockImageFilter< TImage>
::SetNumberOfPartitions(std::vector<unsigned int> input)
{
  this->numOfPartitions = input;
}


template< class TImage>
void BlockImageFilter< TImage>
::GenerateData()
{
  typename TImage::ConstPointer input = this->GetInput();
  typename TImage::Pointer output = this->GetOutput();

  //Allocation
  this->AllocateOutputs();
  itk::ImageAlgorithm::Copy(input.GetPointer(), output.GetPointer(), output->GetRequestedRegion(),
                       output->GetRequestedRegion() );
  output->SetOrigin( input->GetOrigin() );
  output->SetSpacing( input->GetSpacing() );
  output->SetDirection( input->GetDirection() );

  unsigned int Dimension = input->ImageDimension;
  unsigned int NBlocks = 1;
  std::vector<unsigned int> blkSize(Dimension), cumSize(Dimension);
  cumSize[0] = 1;
  typename TImage::SizeType imgSize = input->GetLargestPossibleRegion().GetSize();
  for (unsigned int d = 0; d<Dimension; d++ )
  {
    blkSize[d] = imgSize[d]/this->numOfPartitions[d];
    NBlocks = NBlocks*this->numOfPartitions[d];
    if(d>0)
      cumSize[d] = cumSize[d-1]*this->numOfPartitions[d-1];
  }

  //Find image range:
  typedef itk::MinimumMaximumImageCalculator< TImage > MinMaxCalculatorType;
  typename MinMaxCalculatorType::Pointer minMaxCalculator = MinMaxCalculatorType::New();
  minMaxCalculator->SetImage(input);
  minMaxCalculator->Compute();
  typename TImage::PixelType imgMin = minMaxCalculator->GetMinimum();
  typename TImage::PixelType imgMax = minMaxCalculator->GetMaximum();
  float range = static_cast< float >( imgMax-imgMin );

  itk::ImageRegionIteratorWithIndex< TImage > it(output,output->GetLargestPossibleRegion());
  itk::ImageRegionConstIterator< TImage > imgIt(input,input->GetLargestPossibleRegion());
  for( imgIt.GoToBegin(), it.GoToBegin(); !it.IsAtEnd(); ++it, ++imgIt)
  {
    typename TImage::IndexType idx = it.GetIndex();
    unsigned int pos = 0;
    for(unsigned int d=0; d<Dimension; d++)
    {
      idx[d] = idx[d]/blkSize[d];
      idx[d] = (idx[d]>this->numOfPartitions[d]-1)?(this->numOfPartitions[d]-1):idx[d];
      pos += idx[d]*cumSize[d];
    }
    it.Set( float(pos)*range + imgIt.Get() );
  }







  /*
  itk::Index<2> cornerPixel = input->GetLargestPossibleRegion().GetIndex();
  typename TImage::PixelType newValue = 3;

  this->AllocateOutputs();
  itk::ImageAlgorithm::Copy(input.GetPointer(), output.GetPointer(), output->GetRequestedRegion(),
                       output->GetRequestedRegion() );

		       output->SetPixel( cornerPixel, newValue );*/
}

}// end namespace


#endif
