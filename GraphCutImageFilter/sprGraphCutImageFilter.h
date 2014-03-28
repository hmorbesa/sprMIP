#ifndef __sprGraphCutImageFilter_h
#define __sprGraphCutImageFilter_h

#include "itkImageToImageFilter.h"
#include "sprBoostGraph.h"

namespace spr
{
template< typename TImage, typename TMask>
class GraphCutImageFilter : public itk::ImageToImageFilter< TImage, TImage >
{
public:
  /** Standard class typedefs. */
  typedef GraphCutImageFilter             Self;
  typedef itk::ImageToImageFilter< TImage, TImage > Superclass;
  typedef itk::SmartPointer< Self >        Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GraphCutImageFilter, ImageToImageFilter);

  /** The image to be inpainted in regions where the mask is white.*/
  void SetInputImage(const TImage* image);

  /** The mask to be inpainted. White pixels will be inpainted, black pixels will be passed through to the output.*/
  void SetInputMask(const TMask* mask);

  /** The maximum number of iterations allowed**/
  void SetMaximumNumberOfIterations (unsigned int _arg)
  {
    this->MaximumNumberOfIterations = _arg;
  }

  void SetValueTolerance( double _arg)
  {
    this->ValueTolerance = _arg;
  }

  void SetVerbose( bool _arg=true )
  {
    this->verbose = _arg;
  }

  void EnhanceImage( bool _arg=true )
  {
    this->enhance = _arg;
  }


protected:
  GraphCutImageFilter();
  ~GraphCutImageFilter(){}

  typename TImage::ConstPointer GetInputImage();
  typename TMask::ConstPointer GetInputMask();
  
  /** Does the real work. */
  virtual void GenerateData();

private:
  GraphCutImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
  unsigned int MaximumNumberOfIterations;
  double ValueTolerance;
  bool verbose, enhance;

};
} //namespace SPR


#ifndef ITK_MANUAL_INSTANTIATION
#include "sprGraphCutImageFilter.txx"
#endif


#endif // __sprGraphCutImageFilter_h
