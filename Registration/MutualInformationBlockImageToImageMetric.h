#ifndef __MutualInformationBlockImageToImageMetric_h
#define __MutualInformationBlockImageToImageMetric_h

#include "itkImageToImageMetric.h"
#include "itkPoint.h"
#include "itkIndex.h"
#include "itkBSplineDerivativeKernelFunction.h"
#include "itkArray2D.h"

#include "itkSimpleFastMutexLock.h"

namespace itk
{

  template <typename TFixedImage, typename TMovingImage>

  class MutualInformationBlockImageToImageMetric:
  public ImageToImageMetric<TFixedImage, TMovingImage>
{
public:

	  typedef MutualInformationBlockImageToImageMetric                         Self;
	  typedef ImageToImageMetric<TFixedImage, TMovingImage>                  Superclass;
                               
    
	typedef SmartPointer<Self>                                               Pointer;
    typedef SmartPointer<const Self>                                         ConstPointer;

	itkNewMacro(Self);
	itkTypeMacro(MattesMutualInformationImageToImageMetric,ImageToImageMetric);






	//** Superclass types */
	
	 typedef typename Superclass::TransformType                  TransformType;
     typedef typename Superclass::TransformPointer               TransformPointer;
     typedef typename Superclass::TransformJacobianType          TransformJacobianType;
     typedef typename Superclass::InterpolatorType               InterpolatorType;
     typedef typename Superclass::MeasureType                    MeasureType;
     typedef typename Superclass::DerivativeType                 DerivativeType;
     typedef typename Superclass::ParametersType                 ParametersType;
     typedef typename Superclass::FixedImageType                 FixedImageType;
     typedef typename Superclass::MovingImageType                MovingImageType;
     typedef typename Superclass::MovingImagePointType           MovingImagePointType;
     typedef typename Superclass::FixedImageConstPointer         FixedImageConstPointer;
     typedef typename Superclass::MovingImageConstPointer        MovingImageConstPointer;
     typedef typename Superclass::BSplineTransformWeightsType    BSplineTransformWeightsType;
     typedef typename Superclass::BSplineTransformIndexArrayType BSplineTransformIndexArrayType;

     typedef typename Superclass::CoordinateRepresentationType CoordinateRepresentationType;
     typedef typename Superclass::FixedImageSampleContainer    FixedImageSampleContainer;
     typedef typename Superclass::ImageDerivativesType         ImageDerivativesType;
     typedef typename Superclass::WeightsValueType             WeightsValueType;
     typedef typename Superclass::IndexValueType               IndexValueType;

     typedef typename FixedImageType::OffsetValueType OffsetValueType;


	 itkStaticConstMacro(MovingImageDimension, unsigned int,MovingImageType::ImageDimension);
	 
	 virtual void Initialize(void) 
		 throw ( ExceptionObject );

	 

	 MeasureType GetValue(const ParametersType & parameters) const;

     void GetDerivative(const ParametersType & parameters, DerivativeType & Derivative) const;
     
	 void GetValueAndDerivative(const ParametersType & parameters, MeasureType & Value, DerivativeType & Derivative) const;

	 itkSetClampMacro( NumberOfHistogramBins, SizeValueType,
                    5, NumericTraits<SizeValueType>::max() );
     itkGetConstReferenceMacro(NumberOfHistogramBins, SizeValueType);


	// itkSetMacro(UseExplicitPDFDerivatives, bool);
    // itkGetConstReferenceMacro(UseExplicitPDFDerivatives, bool);
    // itkBooleanMacro(UseExplicitPDFDerivatives);

	  typedef double PDFValueType; //NOTE:  floating point precision is not as stable.  Double precision proves faster and more robust in real-world testing.

      /** Typedef for the joint PDF and PDF derivatives are stored as ITK Images. */
      typedef Image<PDFValueType, 2> JointPDFType;
      typedef Image<PDFValueType, 3> JointPDFDerivativesType;


	  void SetNumberOfHistogramBins(int nHB);
    

protected: 

	MutualInformationBlockImageToImageMetric();
	 
	 virtual ~MutualInformationBlockImageToImageMetric();
	 void PrintSelf(std::ostream & os, Indent indent) const;
	
	

private:

	 
	 void operator=(const Self &);

	 typedef JointPDFType::IndexType             JointPDFIndexType;
     typedef JointPDFType::PixelType             JointPDFValueType;
     typedef JointPDFType::RegionType            JointPDFRegionType;
     typedef JointPDFType::SizeType              JointPDFSizeType;
     typedef JointPDFDerivativesType::IndexType  JointPDFDerivativesIndexType;
     typedef JointPDFDerivativesType::PixelType  JointPDFDerivativesValueType;
     typedef JointPDFDerivativesType::RegionType JointPDFDerivativesRegionType;
     typedef JointPDFDerivativesType::SizeType   JointPDFDerivativesSizeType;

	// Revisar
	   
	 /** Variables to define the marginal and joint histograms. */
      SizeValueType m_NumberOfHistogramBins;
	  int           m_numOfPart;
	  
      PDFValueType  m_MovingImageNormalizedMin;
      PDFValueType  m_FixedImageNormalizedMin;
      PDFValueType  m_FixedImageTrueMin;
      PDFValueType  m_FixedImageTrueMax;
      PDFValueType  m_MovingImageTrueMin;
      PDFValueType  m_MovingImageTrueMax;
      PDFValueType  m_FixedImageBinSize;
      PDFValueType  m_MovingImageBinSize;



     //virtual void ComputeResults(void) const;

  };
}

#ifndef ITK_MANUAL_INSTANTIATION
#include "MutualInformationBlockImageToImageMetric.hxx"
#endif

#endif
