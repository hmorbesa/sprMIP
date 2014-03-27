/*=========================================================================
 *
 * Most of this code was taken from the example for itkLabelOverlapMeasuresImageFilter
 * provided in http://www.insight-journal.org/browse/publication/707
 *
 *=========================================================================*/
#include <iomanip>

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkAddImageFilter.h"
#include "itkLabelOverlapMeasuresImageFilter.h"

const    unsigned int   Dimension = 3;
typedef  unsigned short LabelSpelType;

//  The types of the input images are instantiated by the following lines.
//
typedef itk::Image< LabelSpelType, Dimension >  LabelImageType;
typedef itk::ImageFileReader< LabelImageType >  LabelImageReaderType;
typedef itk::AddImageFilter< LabelImageType >   AddImageType;
typedef itk::LabelOverlapMeasuresImageFilter<
            LabelImageType >                    LabelMeasuresType;

int main( int argc, char *argv[] )
{
  if( argc < 3 )
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage:\n" << argv[0];
    std::cerr << " sourceImg targetImg [use_background]\n";
    std::cerr << "Optional:\n"
              << "\tuse_background [bool]: Compute measures for background also";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  LabelImageReaderType::Pointer srcReader   = LabelImageReaderType::New();
  LabelImageReaderType::Pointer tgtReader   = LabelImageReaderType::New();

  srcReader->SetFileName(argv[1]);
  tgtReader->SetFileName(argv[2]);

  try
  {
    srcReader->Update();
    tgtReader->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Exception thrown while reading the images" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  LabelMeasuresType::Pointer performance = LabelMeasuresType::New();

  if( argc > 3 )
  {
    AddImageType::Pointer srcAdd = AddImageType::New();
    AddImageType::Pointer tgtAdd = AddImageType::New();
    srcAdd->SetInput1( srcReader->GetOutput() );
    srcAdd->SetConstant2( 1 );
    srcAdd->Update();
    tgtAdd->SetInput1( tgtReader->GetOutput() );
    tgtAdd->SetConstant2( 1 );
    tgtAdd->Update();
    performance->SetSourceImage( srcAdd->GetOutput() );
    performance->SetTargetImage( tgtAdd->GetOutput() );
  }
  else
  {
    performance->SetSourceImage( srcReader->GetOutput() );
    performance->SetTargetImage( tgtReader->GetOutput() );
  }

  performance->Update();

  /*std::cout << "                                          "
            << "************ All Labels *************" << std::endl;
            */
  std::cout << std::setw( 10 ) << "Label"
    << std::setw( 17 ) << "Total"
    << std::setw( 17 ) << "Union (jaccard)"
    << std::setw( 17 ) << "Mean (dice)"
    << std::setw( 17 ) << "Volume sim."
    << std::setw( 17 ) << "False negative"
    << std::setw( 17 ) << "False positive" << std::endl;

  std::cout << std::setw( 10 ) << " ";
  std::cout << std::setw( 17 ) << performance->GetTotalOverlap();
  std::cout << std::setw( 17 ) << performance->GetUnionOverlap();
  std::cout << std::setw( 17 ) << performance->GetMeanOverlap();
  std::cout << std::setw( 17 ) << performance->GetVolumeSimilarity();
  std::cout << std::setw( 17 ) << performance->GetFalseNegativeError();
  std::cout << std::setw( 17 ) << performance->GetFalsePositiveError();
  std::cout << std::endl;

  /*std::cout << "                                       "
            << "************ Individual Labels *************" << std::endl;
  std::cout << std::setw( 10 ) << "Label"
            << std::setw( 17 ) << "Target"
            << std::setw( 17 ) << "Union (jaccard)"
            << std::setw( 17 ) << "Mean (dice)"
            << std::setw( 17 ) << "Volume sim."
            << std::setw( 17 ) << "False negative"
            << std::setw( 17 ) << "False positive" << std::endl;
            */

  LabelMeasuresType::MapType labelMap = performance->GetLabelSetMeasures();
  LabelMeasuresType::MapType::const_iterator it;
  for( it = labelMap.begin(); it != labelMap.end(); ++it )
    {
    if( (*it).first == 0 )
      {
      continue;
      }

    int label = (*it).first;

    std::cout << std::setw( 10 ) << label;
    std::cout << std::setw( 17 ) << performance->GetTargetOverlap( label );
    std::cout << std::setw( 17 ) << performance->GetUnionOverlap( label );
    std::cout << std::setw( 17 ) << performance->GetMeanOverlap( label );
    std::cout << std::setw( 17 ) << performance->GetVolumeSimilarity( label );
    std::cout << std::setw( 17 ) << performance->GetFalseNegativeError( label );
    std::cout << std::setw( 17 ) << performance->GetFalsePositiveError( label );
    std::cout << std::endl;
    }

  return EXIT_SUCCESS;
}
