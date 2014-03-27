#ifndef __sprMRnormalization_h
#define __sprMRnormalization_h

#include "itkImageToImageFilter.h"

namespace spr
{
template< class TImage>
class mrNormalization:public itk::ImageToImageFilter< TImage, TImage >
{
public:

  /** Standard class typedefs. */
  typedef mrNormalization             Self;
  typedef itk::ImageToImageFilter< TImage, TImage > Superclass;
  typedef itk::SmartPointer< Self >        Pointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(mrNormalization, ImageToImageFilter);

  itkSetMacro( Variable, double );
  itkGetMacro( Variable, double);

  //Method for updating the lowest allowed quantile
  void SetLowestQuantile(double value)
  {
    this->tl = value;
  }

  //Method for updating the highest allowed quantile
  void SetHighestQuantile(double value)
  {
    this->tu = value;
  }

protected:
  mrNormalization(){this->tl=0.02; this->tu=0.98;}
  ~mrNormalization(){}

  double m_Variable;

  /** Does the real work. */
  virtual void GenerateData();

private:
  mrNormalization(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
  double tl;                     //lowest quantile
  double tu;                     //highest quantile

};
} //namespace SPR


#ifndef ITK_MANUAL_INSTANTIATION
#include "sprMRnormalization.hxx"
#endif

#endif // __sprMRnormalization_h
