/*=========================================================================
 *
 *
 *
 *
 *=========================================================================*/

#include "vector"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkVectorImage.h"

const    unsigned int   Dimension = 3;
typedef  unsigned int   LabelSpelType;


//  The types of the input images are instantiated by the following lines.
//
typedef itk::Image< LabelSpelType, Dimension >          LabelImageType;
typedef itk::ImageFileReader< LabelImageType  >         LabelImageReaderType;
typedef itk::ImageFileWriter< LabelImageType  >         LabelImageWriterType;
typedef std::vector<LabelSpelType>                      VectorContainerType;
typedef VectorContainerType::iterator                   VectorIteratorType;
typedef itk::VectorImage<VectorContainerType,
            Dimension>                                  VotingLabelType;
typedef itk::ImageRegionConstIterator< LabelImageType > ConstLabelIteratorType;
typedef itk::ImageRegionIterator< VotingLabelType >     VotingLabelIteratorType;

void CreateVotingImage(VotingLabelType::Pointer p, LabelImageType::Pointer refImg)
{
  p->SetRegions(refImg->GetLargestPossibleRegion());
  p->SetOrigin(refImg->GetOrigin());
  p->SetSpacing(refImg->GetSpacing());
  p->SetDirection(refImg->GetDirection());
  p->SetVectorLength(2);
  p->Allocate();
}

int main( int argc, char *argv[] )
{
  if( argc < 2 )
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage:\n" << argv[0];
    std::cerr << " listOfImages";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  /****Begin Variable declaration****/
  unsigned int N = 0;                                                       //Number of images in the dataset
  VotingLabelType::Pointer poll           = VotingLabelType::New();         //Image for storing votes
  LabelImageReaderType::Pointer lblReader = LabelImageReaderType::New();
  /****End Variable declaration****/

  //Auxiliar variables:
  LabelImageType::Pointer input;
  std::ifstream file(argv[1]);
  std::string lblName;
  bool configAtlas = true;

  while (getline(file, lblName, '\n'))
  {
    ++N;
    lblReader->SetFileName( lblName );
    try
    {
      lblReader->Update();
    }
    catch( itk::ExceptionObject & excp )
    {
      std::cerr << "Exception thrown while reading the image" << std::endl;
      std::cerr << excp << std::endl;
      return EXIT_FAILURE;
    }
    input = lblReader->GetOutput();
    input->DisconnectPipeline();

    if(configAtlas)
    {
      CreateVotingImage(poll, lblReader->GetOutput());
      configAtlas = false;
    }

    //Setup Iterators:
    ConstLabelIteratorType  itLbl(input,input->GetLargestPossibleRegion());
    VotingLabelIteratorType itPoll(poll,poll->GetLargestPossibleRegion());

    //MAIN LOOP:
    for( itPoll.GoToBegin(), itLbl.GoToBegin(); !itLbl.IsAtEnd(); ++itPoll, ++itLbl )
    {
      LabelImageType::InternalPixelType c                 = itLbl.Get();
      itk::VariableLengthVector<VectorContainerType> vote = itPoll.Get();

      VectorIteratorType it = std::find(vote[0].begin(), vote[0].end(), c);

      if(it==vote[0].end()) //If there's not that label at that spel
      {
        std::cout << vote.GetSize() << std::endl;
        vote[0].push_back(c);//Include it...
        vote[1];//.push_back(1);//with 1 vote
        std::cout << "5\n";
      }
      else //just add 1 vote
      {
        std::cout << "6\n";
        size_t index = std::distance(vote[0].begin(), it);
        std::cout << "7\n";
        vote[1].at(index) ++ ;//with 1 vote
        std::cout << "8\n";
      }
    }
  }



/*

  ImageType::IndexType start;
  start.Fill(0);
  ImageType::SizeType size;
  size.Fill(2);
  ImageType::RegionType region(start,size);
  ImageType::Pointer image = ImageType::New();
  image->SetRegions(region);
  image->SetVectorLength(2);
  image->Allocate();

  ImageType::IndexType pixelIndex;
  ImageType::PixelType pixelValue;
  pixelIndex[0] = 1;
  pixelIndex[1] = 1;
  pixelIndex[2] = 1;
  pixelValue = image->GetPixel(pixelIndex);
  pixelValue[0].push_back(0);
  pixelValue[0].push_back(2);
  pixelValue[1].push_back(1);
  pixelValue[1].push_back(3);
  image->SetPixel(pixelIndex, pixelValue);

  for ( int i = 0; i < pixelValue[0].size(); i++)
  {
    std::cout << pixelValue[0][i] << " " << pixelValue[1][i] << std::endl;
  }

  /*
  ImageType::IndexType pixelIndex;
  pixelIndex[0] = 1;
  pixelIndex[1] = 1;

  ImageType::PixelType pixelValue = image->GetPixel(pixelIndex);
  std::cout << "pixel (1,1) = " << pixelValue << std::endl;

  typedef itk::VariableLengthVector<double> VariableVectorType;
  VariableVectorType variableLengthVector;
  variableLengthVector.SetSize(2);
  variableLengthVector[0] = 1.1;
  variableLengthVector[1] = 2.2;

  image->SetPixel(pixelIndex, variableLengthVector);

  std::cout << "pixel (1,1) = " << pixelValue << std::endl;


*/





  return EXIT_SUCCESS;
}
