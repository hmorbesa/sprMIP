/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "sprMRnormalization.h"

#include "itkAddImageFilter.h"
#include "itkDivideImageFilter.h"
#include "itkBinaryThresholdImageFilter.h"
#include "itkVectorImageToImageAdaptor.h"

// The types of each one of the components in the registration methods should
// be instantiated. First, we select the image dimension and the type for
// representing image pixels.
//
const    unsigned int    Dimension = 3;
typedef  float           IntensitySpelType;
typedef  float           LabelSpelType;
const double PI = std::atan(1.0)*4;

//  The types of the input images are instantiated by the following lines.
//
typedef itk::Image< IntensitySpelType, Dimension >          IntensityImageType;
typedef itk::Image< LabelSpelType, Dimension >              LabelImageType;
typedef itk::VectorImage< IntensitySpelType, Dimension >    AtlasType;

typedef itk::ImageFileReader< IntensityImageType  >         IntensityImageReaderType;
typedef itk::ImageFileReader< AtlasType  >                  AtlasReaderType;
typedef itk::ImageFileReader< LabelImageType  >             LabelImageReaderType;
typedef itk::ImageFileWriter< AtlasType  >                  AtlasWriterType;
typedef itk::ImageFileWriter< LabelImageType  >             LabelImageWriterType;

typedef spr::mrNormalization<IntensityImageType>            MRNormalizationFilterType;

typedef itk::ImageRegionConstIterator< IntensityImageType > ConstIntensityIteratorType;
typedef itk::ImageRegionConstIterator< AtlasType >          ConstAtlasIteratorType;
//typedef itk::ImageRegionIterator< IntensityImageType >      IntensityIteratorType;
typedef itk::ImageRegionIterator< LabelImageType >          LabelIteratorType;

typedef itk::ImageRegionConstIteratorWithIndex< AtlasType > ConstAtlasIteratorIndexType;


void CreateLabelImage(LabelImageType::Pointer p, IntensityImageType::Pointer refImg)
{
  p->SetRegions(refImg->GetLargestPossibleRegion());
  p->SetOrigin(refImg->GetOrigin());
  p->SetSpacing(refImg->GetSpacing());
  p->SetDirection(refImg->GetDirection());
  p->Allocate();
  p->FillBuffer(static_cast<LabelImageType::PixelType>(0));
}


int main( int argc, char *argv[] )
{
  if( argc < 6 )
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage:\n" << argv[0];
    std::cerr << " priors mu sigma queryImg outImg";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  AtlasReaderType::Pointer prReader             = AtlasReaderType::New();
  prReader->SetFileName(argv[1]);

  AtlasReaderType::Pointer muReader             = AtlasReaderType::New();
  muReader->SetFileName(argv[2]);

  AtlasReaderType::Pointer sigmaReader          = AtlasReaderType::New();
  sigmaReader->SetFileName(argv[3]);

  IntensityImageReaderType::Pointer imgReader   = IntensityImageReaderType::New();
  imgReader->SetFileName(argv[4]);

  try
  {
    prReader->Update();
    muReader->Update();
    sigmaReader->Update();
    imgReader->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Exception thrown while reading the images" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  unsigned int C = prReader->GetOutput()->GetNumberOfComponentsPerPixel();

  //Intensity normalization
  MRNormalizationFilterType::Pointer mrNormalization = MRNormalizationFilterType::New();
  mrNormalization->SetInput( imgReader->GetOutput() );
  mrNormalization->Update();

  LabelImageType::Pointer lblImage = LabelImageType::New();
  CreateLabelImage(lblImage,mrNormalization->GetOutput());

  //Setup Iterators:
  ConstIntensityIteratorType  itImg(mrNormalization->GetOutput(), mrNormalization->GetOutput()->GetLargestPossibleRegion());
  ConstAtlasIteratorType      itPr(prReader->GetOutput(), prReader->GetOutput()->GetLargestPossibleRegion());
  ConstAtlasIteratorType      itMu(muReader->GetOutput(),muReader->GetOutput()->GetLargestPossibleRegion());
  ConstAtlasIteratorType      itSig(sigmaReader->GetOutput(),sigmaReader->GetOutput()->GetLargestPossibleRegion());
  LabelIteratorType           itLbl(lblImage,lblImage->GetLargestPossibleRegion());

  //Main loop:
  for( itSig.GoToBegin(), itMu.GoToBegin(),
         itPr.GoToBegin(), itLbl.GoToBegin(), itImg.GoToBegin();
       !itImg.IsAtEnd();
       ++itSig, ++itMu, ++itPr, ++itLbl, ++itImg )
  {
    IntensityImageType::InternalPixelType x  = itImg.Get();
    LabelImageType::InternalPixelType     l;
    AtlasType::PixelType                  pr = itPr.Get();
    AtlasType::PixelType                  mu = itMu.Get();
    AtlasType::PixelType                  sig = itSig.Get();
    IntensityImageType::InternalPixelType v, v_ant = -1;

    for (unsigned int c=0; c<C; c++)
    {
      sig[c] += 1e-6;
      v = std::exp(-0.5*std::pow(x-mu[c],2.0)/sig[c])/std::sqrt(2*PI*sig[c])*pr[c];
      if (v>v_ant)
      {
        v_ant = v;
        l = c;
      }
    }
    itLbl.Set(l);
  }

  LabelImageWriterType::Pointer lblWriter = LabelImageWriterType::New();
  lblWriter->SetInput(lblImage);
  lblWriter->SetFileName(argv[5]);

  try
  {
    lblWriter->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Exception thrown while writing the image" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Done!\n";

  return EXIT_SUCCESS;
}
