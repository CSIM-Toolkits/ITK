#ifndef __itkIsotropicAnomalousDiffusionImageFilter_h
#define __itkIsotropicAnomalousDiffusionImageFilter_h
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkNumericTraits.h"

namespace itk
{

template< typename TInputImage, typename TOutputImage >
class ITK_EXPORT itkIsotropicAnomalousDiffusionImageFilter:
  public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
  /** Extract dimension from input and output image. */
  itkStaticConstMacro(InputImageDimension, unsigned int,
                      TInputImage::ImageDimension);
  itkStaticConstMacro(OutputImageDimension, unsigned int,
                      TOutputImage::ImageDimension);

  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage  InputImageType;
  typedef TOutputImage OutputImageType;

  /** Standard class typedefs. */
  typedef itkIsotropicAnomalousDiffusionImageFilter                             Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType >                 Superclass;
  typedef SmartPointer< Self >                                                  Pointer;
  typedef SmartPointer< const Self >                                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods). */
  itkTypeMacro(itkIsotropicAnomalousDiffusionImageFilter, ImageToImageFilter)

  typedef typename InputImageType::PixelType                 InputPixelType;
  typedef typename OutputImageType::PixelType                OutputPixelType;
  typedef typename NumericTraits< InputPixelType >::RealType InputRealType;

  itkSetMacro(GeneralizedDiffusion, double)
  itkSetMacro(Iterations, int)
  itkSetMacro(TimeStep, double)
  itkSetMacro(Q, double)
  itkGetMacro(GeneralizedDiffusion, double)
  itkGetMacro(Iterations, int)
  itkGetMacro(TimeStep, double)
  itkGetMacro(Q, double)

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro( InputHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< InputPixelType > ) );
  itkConceptMacro( SameDimensionCheck,
                   ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );

#endif

protected:
  itkIsotropicAnomalousDiffusionImageFilter();
  virtual ~itkIsotropicAnomalousDiffusionImageFilter() {}
    double m_GeneralizedDiffusion;
    int m_Iterations;
    double m_TimeStep;
    double m_Q;
  void GenerateData();

private:
  itkIsotropicAnomalousDiffusionImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
  void TimeStepTestStability();
  double GeneralizedDiffCurve();
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkIsotropicAnomalousDiffusionImageFilter.hxx"
#endif

#endif
