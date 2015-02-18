#ifndef __sprBlockImageFilter_h
#define __sprBlockImageFilter_h

#include "itkImageToImageFilter.h"

namespace spr
{
template< class TImage>
  class BlockImageFilter:public itk::ImageToImageFilter< TImage, TImage >
{
public:

  /** Standard class typedefs. */
  typedef BlockImageFilter             Self;
  typedef itk::ImageToImageFilter< TImage, TImage > Superclass;
  typedef itk::SmartPointer< Self >        Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(BlockImageFilter, ImageToImageFilter);

  itkSetMacro( Variable, double );
  itkGetMacro( Variable, double);

  void SetNumberOfPartitions(std::vector<unsigned int> input);

protected:
  BlockImageFilter(){}
  ~BlockImageFilter(){}

  /** Does the real work. */
  virtual void GenerateData();

  double m_Variable;
  std::vector<unsigned int> numOfPartitions;

private:
  BlockImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented

};
} //namespace SPR


#ifndef ITK_MANUAL_INSTANTIATION
#include "sprBlockImageFilter.hxx"
#endif


#endif // __sprBlockImageFilter_h
