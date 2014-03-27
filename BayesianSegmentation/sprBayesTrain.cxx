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

//  The types of the input images are instantiated by the following lines.
//
typedef itk::Image< IntensitySpelType, Dimension >  IntensityImageType;
typedef itk::Image< LabelSpelType, Dimension >  LabelImageType;
typedef itk::VectorImage< IntensitySpelType, Dimension > AtlasType;

typedef itk::ImageFileReader< IntensityImageType  > IntensityImageReaderType;
typedef itk::ImageFileReader< LabelImageType  > LabelImageReaderType;
typedef itk::ImageFileWriter< IntensityImageType  > IntensityImageWriterType;

typedef itk::AddImageFilter<IntensityImageType> AddIntensityFilterType;

typedef spr::mrNormalization<IntensityImageType> MRNormalizationFilterType;

typedef itk::ImageRegionConstIterator< IntensityImageType >   ConstIntensityIteratorType;
typedef itk::ImageRegionIterator< IntensityImageType >        IntensityIteratorType;
typedef itk::ImageRegionConstIterator< LabelImageType >       ConstLabelIteratorType;
typedef itk::ImageRegionIterator< LabelImageType >            LabelIteratorType;
typedef itk::ImageRegionIterator< AtlasType >                 AtlasIteratorType;


void CreateAtlas(AtlasType::Pointer p, IntensityImageType::Pointer refImg, unsigned int numClass)
{
  p->SetRegions(refImg->GetLargestPossibleRegion());
  p->SetOrigin(refImg->GetOrigin());
  p->SetSpacing(refImg->GetSpacing());
  p->SetDirection(refImg->GetDirection());
  p->SetVectorLength(numClass);
  p->Allocate();

  itk::VariableLengthVector< IntensitySpelType > v0;
  v0.SetSize(numClass);
  for (unsigned int c=0; c<numClass; c++)
    v0[c] = static_cast<AtlasType::PixelType::ValueType>(0.0);

  p->FillBuffer(v0);
}

void CreateImage(IntensityImageType::Pointer p, IntensityImageType::Pointer refImg)
{
  p->SetRegions(refImg->GetLargestPossibleRegion());
  p->SetOrigin(refImg->GetOrigin());
  p->SetSpacing(refImg->GetSpacing());
  p->SetDirection(refImg->GetDirection());
  p->Allocate();
  p->FillBuffer(static_cast<IntensityImageType::PixelType>(0.0));
}

int main( int argc, char *argv[] )
{
  if( argc < 4 )
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage:\n" << argv[0];
    std::cerr << " listOfLabels listOfIntensityImages listOfLabelImages";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  /****Begin Variable declaration****/
  IntensityImageType::PixelType N = 0;                              //Number of images in the dataset
  unsigned int C;                                                   //Number of considered classes
  std::vector<LabelImageType::PixelType> L;                         //Label set: size(L) = C
  AtlasType::Pointer mu = AtlasType::New();                         //Conditional: P(x|c) = exp(-(x-mu_c)^2/s_c^2)
  AtlasType::Pointer m2 = AtlasType::New();                         //Temporal variable
  AtlasType::Pointer Pr = AtlasType::New();                         //Priors: P(c)
  IntensityImageType::Pointer tpt = IntensityImageType::New();      // Template image: intensity dataset average
  IntensityImageReaderType::Pointer imgReader = IntensityImageReaderType::New();
  LabelImageReaderType::Pointer lblReader = LabelImageReaderType::New();
  /****End Variable declaration****/

  //Loading the labels
  std::ifstream file1(argv[1]);
  std::string stuff;
  while (getline(file1, stuff, '\n'))
  {
    L.push_back(atoi(stuff.c_str()));
  }
  file1.close();

  C = L.size();

  //Some configurations:  
  MRNormalizationFilterType::Pointer mrNormalization = MRNormalizationFilterType::New();
  mrNormalization->SetInput( imgReader->GetOutput() );

  IntensityImageType::Pointer input;

  //Auxiliar variables:
  std::ifstream file2(argv[2]);
  std::ifstream file3(argv[3]);
  std::string imgName;
  std::string lblName;
  bool configAtlas = true;

  while (getline(file2, imgName, '\n') && getline(file3, lblName, '\n'))
  {
    ++N;
    //std::cout << ++N << std::endl;

    imgReader->SetFileName( imgName );
    lblReader->SetFileName( lblName );

    try
    {
      imgReader->Update();
      lblReader->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Exception thrown while reading the image" << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }

    mrNormalization->Update();
    input = mrNormalization->GetOutput();
    input->DisconnectPipeline();

    if(configAtlas)
    {
      CreateAtlas(Pr, imgReader->GetOutput(), C);      
      CreateAtlas(m2, imgReader->GetOutput(), C);
      CreateAtlas(mu, imgReader->GetOutput(), C);
      CreateImage(tpt, imgReader->GetOutput());

      IntensityImageWriterType::Pointer imgWriter1 = IntensityImageWriterType::New();
      imgWriter1->SetInput(imgReader->GetOutput());
      imgWriter1->SetFileName("img.nii");
      imgWriter1->Update();


      configAtlas = false;
    }

    //Setup Iterators:
    ConstIntensityIteratorType  itImg(input,input->GetLargestPossibleRegion());
    IntensityIteratorType       itTpt(tpt,tpt->GetLargestPossibleRegion());
    ConstLabelIteratorType      itLbl(lblReader->GetOutput(),lblReader->GetOutput()->GetLargestPossibleRegion());
    AtlasIteratorType           itPr(Pr,Pr->GetLargestPossibleRegion());
    AtlasIteratorType           itMu(mu,mu->GetLargestPossibleRegion());
    AtlasIteratorType           itM2(m2,m2->GetLargestPossibleRegion());

    for( itM2.GoToBegin(), itMu.GoToBegin(),
           itPr.GoToBegin(), itLbl.GoToBegin(), itTpt.GoToBegin(), itImg.GoToBegin();
         !itImg.IsAtEnd();
         ++itM2, ++itMu, ++itPr, ++itLbl, ++itTpt, ++itImg )
    {
      IntensityImageType::InternalPixelType x   = itImg.Get();
      LabelImageType::InternalPixelType     c   = itLbl.Get();
      AtlasType::PixelType                  sPr = itPr.Get();
      AtlasType::PixelType                  sMu = itMu.Get();      
      AtlasType::PixelType                  sM2 = itM2.Get();

      sPr[c]++;

      double delta = x - sMu[c];
      sMu[c] += delta/sPr[c];
      sM2[c] += delta*(x - sMu[c]);

      itTpt.Set( itTpt.Get() + x );  //Updating the template

      itPr.Set( sPr );               //Updating the priors
      itMu.Set( sMu );               //Updating the means      
      itM2.Set( sM2 );               //Updating temporal variables
    }
  }
  file2.close();
  file3.close();

//Last iteration to finally compute priors and var
  IntensityIteratorType       itTpt(tpt,tpt->GetLargestPossibleRegion());
  AtlasIteratorType           itPr(Pr,Pr->GetLargestPossibleRegion());  
  AtlasIteratorType           itM2(m2,m2->GetLargestPossibleRegion());
  for( itTpt.GoToBegin(), itM2.GoToBegin(), itPr.GoToBegin(); !itPr.IsAtEnd();
       ++itTpt, ++itM2, ++itPr )
  {
    AtlasType::PixelType                  sPr = itPr.Get();    
    AtlasType::PixelType                  sM2 = itM2.Get();
    itTpt.Set( itTpt.Get()/N );
    for(unsigned int c=0; c<C; c++)
    {      
      if (sPr[c]>1)
        sM2[c] = sM2[c]/(sPr[c]-1);
      sPr[c] = sPr[c]/N;
    }
    itPr.Set( sPr );
    itM2.Set( sM2 );
  }

  IntensityImageWriterType::Pointer imgWriter = IntensityImageWriterType::New();
  imgWriter->SetInput(tpt);
  imgWriter->SetFileName("tpt.nii");
  imgWriter->Update();

  itk::ImageFileWriter<AtlasType>::Pointer prWriter = itk::ImageFileWriter<AtlasType>::New();
  prWriter->SetInput(Pr);
  prWriter->SetFileName("Pr.nii");
  prWriter->Update();

  itk::ImageFileWriter<AtlasType>::Pointer muWriter = itk::ImageFileWriter<AtlasType>::New();
  muWriter->SetInput(mu);
  muWriter->SetFileName("mu.nii");
  muWriter->Update();

  itk::ImageFileWriter<AtlasType>::Pointer sigWriter = itk::ImageFileWriter<AtlasType>::New();
  sigWriter->SetInput(m2);
  sigWriter->SetFileName("sigma.nii");
  sigWriter->Update();

  std::cout << "Done!\n";

  return EXIT_SUCCESS;
}
