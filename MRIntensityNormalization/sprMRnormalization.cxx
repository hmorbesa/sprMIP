/*=========================================================================
 * Magnetic Resonance Intensity Normalization as proposed in
 * Smith, S. M. (2002). Fast robust automated brain extraction. Human Brain
 * Mapping, 17(3), 143–55. doi:10.1002/hbm.10062
 *
 * Intensity values are normalized to belong to the interval [0,1] as:
 *             y = min( max( (x-x_l)/(x_u-x_l), 0 ), 1 )
 * where x_l and x_u correspond to robust image intensity minimum and
 * maximum, respectively. x_l is the intensity below which lies 2% of the
 * cumulative histogram. x_u is the intensity above which lies 98% of the
 * cumulative histogram.
 *
 * The calculation of the above extrema allows to ignore small numbers of
 * voxels having widely different values from the rest of the image.
 * Specifically, the robust maximum allows to avoid the influence of a few
 * high intensity “outlier” voxels; for example, the DC spike from image
 * reconstruction, or arteries, which often appear much brighter than the
 * rest of the image.
 *
 * For the sake of generalization, the lowest and higest allowed quantile
 * can be provided as function parameters.
 *
 * David C\'ardenas Pe\~na
 * Signal Processing and Recognition Group
 * Timestamp: 2014-03-19 13:12:09 $
 *=========================================================================*/

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "sprMRnormalization.h"

// Types instantiation:
const    unsigned int    Dimension = 3;
typedef  float           IntensitySpelType;
typedef itk::Image< IntensitySpelType, Dimension >  IntensityImageType;
typedef itk::ImageFileReader< IntensityImageType  > IntensityImageReaderType;
typedef itk::ImageFileWriter< IntensityImageType  > IntensityImageWriterType;

typedef spr::mrNormalization<IntensityImageType> MRNormalizationFilterType;

int main( int argc, char *argv[] )
{
  if( argc < 3 )
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage:\n";
    std::cerr << argv[0] << " input_image output_image\n";
    std::cerr << argv[0] << " input_image output_image lowest_quartile highest_quartile\n";
    return EXIT_FAILURE;
  }

  IntensityImageReaderType::Pointer imgReader = IntensityImageReaderType::New();
  imgReader->SetFileName( argv[1] );
  try
  {
    imgReader->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Exception thrown while reading the image" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  MRNormalizationFilterType::Pointer mrNormalization = MRNormalizationFilterType::New();  
  mrNormalization->SetInput( imgReader->GetOutput() );
  mrNormalization->Update();

  IntensityImageWriterType::Pointer imgWriter = IntensityImageWriterType::New();
  imgWriter->SetFileName( argv[2] );
  imgWriter->SetInput(mrNormalization->GetOutput());

  try
  {
    imgWriter->Update();
  }
  catch( itk::ExceptionObject & excp )
  {
    std::cerr << "Exception thrown while writing the image" << std::endl;
    std::cerr << excp << std::endl;
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
