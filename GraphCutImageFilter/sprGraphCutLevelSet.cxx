
#include "itkWin32Header.h"
#include <iostream>
#include <fstream>
#include "itkNumericTraits.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "sprGraphCutImageFilter.h"

const unsigned int Dimension = 2;

typedef itk::Image< unsigned char, Dimension > ImageType;
typedef itk::ImageFileReader< ImageType >      ReaderType;
typedef itk::ImageFileWriter< ImageType >      WriterType;
typedef spr::GraphCutImageFilter< ImageType,
            ImageType >                        GraphCutSegmenterType;

int main(int argc, char **argv)
{

  //Usage: argv[0] inputImage inputMask outputImage

    //Create and setup the readers
    ReaderType::Pointer imgReader = ReaderType::New();
    imgReader->SetFileName(argv[1]);
    imgReader->Update();
    ReaderType::Pointer mskReader = ReaderType::New();
    mskReader->SetFileName(argv[2]);
    mskReader->Update();

    // Create and setup the graphcut-based segmenter
    GraphCutSegmenterType::Pointer segmenter = GraphCutSegmenterType::New();
    segmenter->SetMaximumNumberOfIterations( atoi(argv[4]) );
    segmenter->SetValueTolerance( atof(argv[5]) );
    segmenter->SetInput( imgReader->GetOutput() );
    segmenter->SetInputMask( mskReader->GetOutput() );
    segmenter->EnhanceImage(false);
    segmenter->SetVerbose();
    segmenter->Update();

    // Create and setup the writer
    WriterType::Pointer writer = WriterType::New();
    writer->SetFileName( argv[3] );
    writer->SetInput(segmenter->GetOutput());
    writer->Update();

  return EXIT_SUCCESS;
}


