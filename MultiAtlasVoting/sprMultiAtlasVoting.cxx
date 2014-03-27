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
typedef itk::ImageRegionIterator< LabelImageType >      LabelIteratorType;
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

void CreateLabelImage(LabelImageType::Pointer p, LabelImageType::Pointer refImg)
{
  p->SetRegions(refImg->GetLargestPossibleRegion());
  p->SetOrigin(refImg->GetOrigin());
  p->SetSpacing(refImg->GetSpacing());
  p->SetDirection(refImg->GetDirection());
  p->Allocate();
//  p->FillBuffer(0);
}

int main( int argc, char *argv[] )
{
  if( argc < 3 )
  {
    std::cerr << "Missing Parameters " << std::endl;
    std::cerr << "Usage:\n" << argv[0];
    std::cerr << " listOfImages outputImage";
    std::cerr << std::endl;
    return EXIT_FAILURE;
  }

  /****Begin Variable declaration****/
  unsigned int N = 0;                                                       //Number of images in the dataset
  VotingLabelType::Pointer poll           = VotingLabelType::New();         //Image for storing votes
  LabelImageReaderType::Pointer lblReader = LabelImageReaderType::New();
  LabelImageType::Pointer output          = LabelImageType::New();          //Output label image
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

    if(configAtlas) //Allocate and configure images
    {
      CreateVotingImage(poll, input);
      CreateLabelImage(output, input);
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
        vote[0].push_back(c);//Include it...
        vote[1].push_back(1);//with 1 vote
      }
      else //just add 1 vote
      {
        size_t index = std::distance(vote[0].begin(), it);
        vote[1].at(index) ++ ;//with 1 vote
      }
    }
  }

  //Setup Iterators:
  LabelIteratorType itLbl(output,output->GetLargestPossibleRegion());
  VotingLabelIteratorType itPoll(poll,poll->GetLargestPossibleRegion());
  std::cout << output->GetLargestPossibleRegion() << std::endl;
  //MAIN LOOP:
  for( itPoll.GoToBegin(), itLbl.GoToBegin(); !itLbl.IsAtEnd(); ++itPoll, ++itLbl )
  {    
    itk::VariableLengthVector<VectorContainerType> vote = itPoll.Get();
    VectorIteratorType it = std::max_element(vote[1].begin(), vote[1].end()); //The most voted label
    size_t index = std::distance(vote[1].begin(), it);
    itLbl.Set(vote[0].at(index));
  }
  //Now save the result
  LabelImageWriterType::Pointer lblWriter = LabelImageWriterType::New();
  lblWriter->SetInput(output);
  lblWriter->SetFileName(argv[2]);
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

  return EXIT_SUCCESS;
}
