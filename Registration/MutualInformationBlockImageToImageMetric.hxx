#ifndef __MutualInformationBlockImageToImageMetric_hxx
#define __MutualInformationBlockImageToImageMetric_hxx

#include "MutualInformationBlockImageToImageMetric.h"
#include "itkMattesMutualInformationImageToImageMetric.h"

#include "itkImageRandomConstIteratorWithIndex.h"
#include "itkImageRegionIterator.h"
#include "itkImageIterator.h"
#include "vnl/vnl_math.h"
#include "itkStatisticsImageFilter.h"

#include "vnl/vnl_vector.h"
#include "vnl/vnl_c_vector.h"

namespace itk
{


template <typename TFixedImage, typename TMovingImage>
MutualInformationBlockImageToImageMetric<TFixedImage, TMovingImage>
::MutualInformationBlockImageToImageMetric()
{
}
 
/*
  m_NumberOfHistogramBins(50),
  m_MovingImageNormalizedMin(0.0),
  m_FixedImageNormalizedMin(0.0),
  m_FixedImageTrueMin(0.0),
  m_FixedImageTrueMax(0.0),
  m_MovingImageTrueMin(0.0),
  m_MovingImageTrueMax(0.0),
  m_FixedImageBinSize(0.0),
  m_MovingImageBinSize(0.0),

  m_CubicBSplineKernel(ITK_NULLPTR),
  m_CubicBSplineDerivativeKernel(ITK_NULLPTR),

  m_PRatioArray(0,0),

  // Initialize memory
  m_MovingImageMarginalPDF(0),

  m_MMIMetricPerThreadVariables(ITK_NULLPTR),

  m_UseExplicitPDFDerivatives(true),
  m_ImplicitDerivativesSecondPass(false)

{

	 this->SetComputeGradient(false); // don't use the default gradient for now
     this->m_WithinThreadPreProcess = true;
     this->m_WithinThreadPostProcess = false;
     this->m_ComputeGradient = false;
}
 */


template <typename TFixedImage, typename TMovingImage>
MutualInformationBlockImageToImageMetric<TFixedImage, TMovingImage>
::~MutualInformationBlockImageToImageMetric()
{

	
}




/*   Functions   */


template <typename TFixedImage, typename TMovingImage>
void MutualInformationBlockImageToImageMetric<TFixedImage, TMovingImage>
	::Initialize(void)
{
    
	this->m_numOfPart = 2;
	


     	 typedef itk::MattesMutualInformationImageToImageMetric < TFixedImage,TMovingImage> MetricType;
	     MetricType::Pointer metric = MetricType::New();
	     metric->SetFixedImage(this->m_FixedImage);
	     metric->SetMovingImage(this->m_MovingImage);
	     metric->SetFixedImageRegion(this->m_FixedImage->GetBufferedRegion() );
	     metric->SetTransform(this->m_Transform);
	     metric->SetNumberOfHistogramBins( this->m_NumberOfHistogramBins );
	     metric->SetInterpolator(this->m_Interpolator);
	     metric->Initialize();


}


template <typename TFixedImage, typename TMovingImage>
void MutualInformationBlockImageToImageMetric<TFixedImage, TMovingImage>
::PrintSelf(std::ostream& os, Indent indent) const
{
	Superclass::PrintSelf(os, indent);
	  
	/*
	//os << indent << "NumberOfHistogramBins: ";
    //os << this->m_NumberOfHistogramBins << std::endl;
	
	os << indent << "FixedImageNormalizedMin: ";
    os << this->m_FixedImageNormalizedMin << std::endl;
    os << indent << "MovingImageNormalizedMin: ";
    os << this->m_MovingImageNormalizedMin << std::endl;
    os << indent << "MovingImageTrueMin: ";
    os << this->m_MovingImageTrueMin << std::endl;
    os << indent << "MovingImageTrueMax: ";
    os << this->m_MovingImageTrueMax << std::endl;
    os << indent << "FixedImageBinSize: ";
    os << this->m_FixedImageBinSize << std::endl;
    os << indent << "MovingImageBinSize: ";
    os << this->m_MovingImageBinSize << std::endl;
    os << indent << "UseExplicitPDFDerivatives: ";
   
	//os << this->m_UseExplicitPDFDerivatives << std::endl;
   // os << indent << "ImplicitDerivativesSecondPass: ";
   // os << this->m_ImplicitDerivativesSecondPass << std::endl;
   
	
  if( this->m_MMIMetricPerThreadVariables && this->m_MMIMetricPerThreadVariables[0].JointPDF.IsNotNull() )
    {
    os << indent << "JointPDF: ";
    os << this->m_MMIMetricPerThreadVariables[0].JointPDF << std::endl;
    }
  if( this->m_MMIMetricPerThreadVariables && this->m_MMIMetricPerThreadVariables[0].JointPDFDerivatives.IsNotNull() )
    {
    os << indent << "JointPDFDerivatives: ";
    os << this->m_MMIMetricPerThreadVariables[0].JointPDFDerivatives;
    }
	*/

}

template <typename TFixedImage, typename TMovingImage>
void MutualInformationBlockImageToImageMetric<TFixedImage, TMovingImage>
 ::SetNumberOfHistogramBins(int nHB)
{
	this->m_NumberOfHistogramBins = nHB;
}






template <typename TFixedImage, typename TMovingImage>
typename MutualInformationBlockImageToImageMetric<TFixedImage, TMovingImage>::MeasureType
 MutualInformationBlockImageToImageMetric<TFixedImage, TMovingImage>
::GetValue(const ParametersType &parameters ) const
{
	
	MeasureType value;
	value = NumericTraits<MeasureType>::Zero;
	/////////////////////////////REGION///////////////////////////////////////
	this->m_numOfPart;
   
	
	for(int i=0;i<this->m_numOfPart;i++)
	{
		for(int j=0;j<this->m_numOfPart;j++)
		{
			for(int k=0;k<this->m_numOfPart;k++)
			{

                 
	                
	               
					//value = (int) 1;// NumericTraits<MeasureType>::Zero;             

                    typedef itk::Image<double,3> ImageType;
	                ImageType::IndexType regionIndex;

                 	ImageType::RegionType region;
	
                	ImageType::RegionType Lregion=this->m_FixedImage->GetBufferedRegion();
	                ImageType::SizeType size=Lregion.GetSize();
	
	                regionIndex[0] = (int) (i*floor(size[0]/this->m_numOfPart));
                    regionIndex[1] = (int) (j*floor(size[1]/this->m_numOfPart));
                    regionIndex[2] = (int) (k*floor(size[2]/this->m_numOfPart));
                    ImageType::SizeType regionSize;
	                regionSize[0] = (int) floor(size[0]/this->m_numOfPart);
                    regionSize[1] = (int) floor(size[1]/this->m_numOfPart);
                    regionSize[2] = (int) floor(size[2]/this->m_numOfPart);   
                    region.SetIndex(regionIndex);
                    region.SetSize(regionSize); 

                    ////////////////////////////////////////////////
 
	                typedef itk::MattesMutualInformationImageToImageMetric < TFixedImage,TMovingImage> MetricType;
	                MetricType::Pointer metric = MetricType::New();
	                metric->SetFixedImage(this->m_FixedImage);
	                metric->SetMovingImage(this->m_MovingImage);
		            //metric->SetFixedImageRegion(this->m_FixedImage->GetBufferedRegion() );
		            metric->SetFixedImageRegion(region);
		            metric->SetTransform(this->m_Transform);
	                metric->SetNumberOfHistogramBins( this->m_NumberOfHistogramBins );
	                metric->SetInterpolator(this->m_Interpolator);
					metric->SetUseAllPixels(true);
	                metric->Initialize();

					value = (metric->GetValue(parameters)+value); //(m_numOfPart*m_numOfPart*m_numOfPart);
					
					//std::cout<<"metric"<< value << metric->GetValue(parameters)  <<std::endl;
                     
	      }
		}

	}
	

    return static_cast<MeasureType>( value/(m_numOfPart*m_numOfPart*m_numOfPart));
  




}


template<typename TFixedImage, typename TMovingImage>
void MutualInformationBlockImageToImageMetric< TFixedImage,TMovingImage>
	::GetDerivative(const ParametersType &parameters, DerivativeType &derivative) const
{
    
	 DerivativeType derivativeacc;

	derivativeacc=DerivativeType(derivative.GetSize());

	/////////////////////////////REGION///////////////////////////////////////
	this->m_numOfPart;
   
	
	for(int i=0;i<this->m_numOfPart;i++)
	{
		for(int j=0;j<this->m_numOfPart;j++)
		{
			for(int k=0;k<this->m_numOfPart;k++)
			{

                 
	                
	               
					//value = (int) 1;// NumericTraits<MeasureType>::Zero;             

                    typedef itk::Image<double,3> ImageType;
	                ImageType::IndexType regionIndex;

                 	ImageType::RegionType region;
	
                	ImageType::RegionType Lregion=this->m_FixedImage->GetBufferedRegion();
	                ImageType::SizeType size=Lregion.GetSize();
	
	                regionIndex[0] = (int)(i*floor(size[0]/this->m_numOfPart));
                    regionIndex[1] = (int)(j*floor(size[1]/this->m_numOfPart));
                    regionIndex[2] = (int)(k*floor(size[2]/this->m_numOfPart));
                    ImageType::SizeType regionSize;
	                regionSize[0] = (int) floor(size[0]/this->m_numOfPart);
                    regionSize[1] = (int) floor(size[1]/this->m_numOfPart);
                    regionSize[2] = (int) floor(size[2]/this->m_numOfPart);   
                    region.SetIndex(regionIndex);
                    region.SetSize(regionSize); 

                    ////////////////////////////////////////////////
 
	                typedef itk::MattesMutualInformationImageToImageMetric < TFixedImage,TMovingImage> MetricType;
	                MetricType::Pointer metric = MetricType::New();
	                metric->SetFixedImage(this->m_FixedImage);
	                metric->SetMovingImage(this->m_MovingImage);
		            //metric->SetFixedImageRegion(this->m_FixedImage->GetBufferedRegion() );
		            metric->SetFixedImageRegion(region);
		            metric->SetTransform(this->m_Transform);
	                metric->SetNumberOfHistogramBins( this->m_NumberOfHistogramBins );
	                metric->SetInterpolator(this->m_Interpolator);
					metric->SetUseAllPixels(true);
	                metric->Initialize();

					
					
					metric->GetDerivative(parameters, derivative);
					
					//std::cout<<"Derivative    "<< derivative[0] << std::endl;

					for(unsigned int parameter=0;  parameter < derivative.GetSize();++parameter)
					{
						derivativeacc[parameter]=derivative[parameter]+derivativeacc[parameter];					
					}
					
					//std::cout<<"Derivativeacc "<< derivativeacc<< std::endl;


					
					/*std::cout<<"derivative"<< derivativeacc<<std::endl;
					std::cout<<"derivative2:"<<derivative<<std::endl;
	                std::cout<<"this->m_NumberOfParameters  " << derivative.GetSize() << std::endl;*/
                     
	      }
		}
      }

                	for(unsigned int parameter;  parameter < derivative.GetSize();++parameter)
					{
						derivative[parameter]=derivativeacc[parameter]/(m_numOfPart*m_numOfPart*m_numOfPart);
					}

}


template<typename TFixedImage, typename TMovingImage>
void MutualInformationBlockImageToImageMetric< TFixedImage,TMovingImage>
::GetValueAndDerivative(const ParametersType & parameters,
                        MeasureType & value,
                        DerivativeType & derivative) const
{
	value = NumericTraits<MeasureType>::Zero;
	
	DerivativeType derivativeacc;
	derivativeacc=DerivativeType(derivative.GetSize());
	
	/////////////////////////////REGION///////////////////////////////////////
	this->m_numOfPart;
   
	
	for(int i=0;i<this->m_numOfPart;i++)
	{
		for(int j=0;j<this->m_numOfPart;j++)
		{
			for(int k=0;k<this->m_numOfPart;k++)
			{

                 
	                
	               
					//value = (int) 1;// NumericTraits<MeasureType>::Zero;             

                    typedef itk::Image<double,3> ImageType;
	                ImageType::IndexType regionIndex;

                 	ImageType::RegionType region;
	
                	ImageType::RegionType Lregion=this->m_FixedImage->GetBufferedRegion();
	                ImageType::SizeType size=Lregion.GetSize();
	
	                regionIndex[0] = (int)(i*floor(size[0]/this->m_numOfPart));
                    regionIndex[1] = (int)(j*floor(size[1]/this->m_numOfPart));
                    regionIndex[2] = (int)(k*floor(size[2]/this->m_numOfPart));
                    ImageType::SizeType regionSize;
	                regionSize[0] = (int) floor(size[0]/this->m_numOfPart);
                    regionSize[1] = (int) floor(size[1]/this->m_numOfPart);
                    regionSize[2] = (int) floor(size[2]/this->m_numOfPart);   
                    region.SetIndex(regionIndex);
                    region.SetSize(regionSize); 

                    ////////////////////////////////////////////////
 
	                typedef itk::MattesMutualInformationImageToImageMetric < TFixedImage,TMovingImage> MetricType;
	                MetricType::Pointer metric = MetricType::New();
	                metric->SetFixedImage(this->m_FixedImage);
	                metric->SetMovingImage(this->m_MovingImage);
		            //metric->SetFixedImageRegion(this->m_FixedImage->GetBufferedRegion() );
		            metric->SetFixedImageRegion(region);
		            metric->SetTransform(this->m_Transform);
	                metric->SetNumberOfHistogramBins( this->m_NumberOfHistogramBins );
	                metric->SetInterpolator(this->m_Interpolator);
					metric->SetUseAllPixels(true);
	                metric->Initialize();

					value = (metric->GetValue(parameters)+value);//(m_numOfPart*m_numOfPart*m_numOfPart);
					

					//std::cout<<"metric"<< value << metric->GetValue(parameters)  <<std::endl;

				//	std::cout<<"metric"<< value<< metric->GetValue(parameters)  <<std::endl;

					metric->GetDerivative(parameters, derivative);

					//std::cout<<"Derivative    "<< derivative[0] << std::endl;

					for(unsigned int parameter=0;  parameter < derivative.GetSize();++parameter)
					{
						derivativeacc[parameter]=derivative[parameter]+derivativeacc[parameter];						
					}					


					//std::cout<<"Derivativeacc "<< derivativeacc<< std::endl;

					//std::cout<<"Derivativeacc "<< derivativeacc[0] << std::endl;

					   /*std::cout<<"derivative"<< derivativeacc<<std::endl;
					   std::cout<<"derivative2:"<<derivative<<std::endl;
	                   std::cout<<"this->m_NumberOfParameters  " << derivative.GetSize()  << std::endl;*/
    
                     
	      }
		}
      }
	value=value/(m_numOfPart*m_numOfPart*m_numOfPart);


                	for(unsigned int parameter=0;  parameter < derivative.GetSize();++parameter)
					{
						derivative[parameter]=derivativeacc[parameter]/(m_numOfPart*m_numOfPart*m_numOfPart);
					}
	
	 
}


};

#endif


