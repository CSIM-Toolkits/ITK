#ifndef __itkAnisotropicAnomalousDiffusionImageFilter_h
#define __itkAnisotropicAnomalousDiffusionImageFilter_h
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkNumericTraits.h"
#include <itkConstNeighborhoodIterator.h>


namespace itk
{

template< typename TInputImage, typename TOutputImage >
class ITK_EXPORT AnisotropicAnomalousDiffusionImageFilter:
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
  typedef AnisotropicAnomalousDiffusionImageFilter                                      Self;
  typedef ImageToImageFilter< InputImageType, OutputImageType >                 Superclass;
  typedef SmartPointer< Self >                                                  Pointer;
  typedef SmartPointer< const Self >                                            ConstPointer;

  typedef typename Superclass::OutputImageRegionType OutputImageRegionType;

  typedef itk::ConstNeighborhoodIterator< OutputImageType >       NeighborIteratorType;

  /** Method for creation through the object factory. */
  itkNewMacro(Self)

  /** Run-time type information (and related methods). */
  itkTypeMacro(AnisotropicAnomalousDiffusionImageFilter, ImageToImageFilter)

  typedef typename InputImageType::PixelType                 InputPixelType;
  typedef typename OutputImageType::PixelType                OutputPixelType;
//  typedef typename NumericTraits< InputPixelType >::RealType InputRealType;

  itkSetMacro(Conductance, float)
  itkSetMacro(Iterations, int)
  itkSetMacro(TimeStep, float)
  itkSetMacro(Q, float)
  itkGetMacro(Conductance, float)
  itkGetMacro(Iterations, int)
  itkGetMacro(TimeStep, float)
  itkGetMacro(Q, float)

#ifdef ITK_USE_CONCEPT_CHECKING
  // Begin concept checking
  itkConceptMacro( InputHasNumericTraitsCheck,
                   ( Concept::HasNumericTraits< InputPixelType > ) );
  itkConceptMacro( SameDimensionCheck,
                   ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );
  itkConceptMacro( InputPixelTypeIsFloatingPointCheck,
                   ( Concept::IsFloatingPoint< InputPixelType > ) );
//  itkConceptMacro( OutputPixelTypeIsFloatingPointCheck,
//                   ( Concept::IsFloatingPoint< OutputPixelType > ) );
  // End concept checking
#endif

protected:
  AnisotropicAnomalousDiffusionImageFilter();
  virtual ~AnisotropicAnomalousDiffusionImageFilter() {}
    float m_GeneralizedDiffusion, m_TimeStep, m_Q, m_Conductance;
    int m_Iterations;
  virtual void ThreadedGenerateData(const OutputImageRegionType &, ThreadIdType);

private:
  AnisotropicAnomalousDiffusionImageFilter(const Self &); //purposely not implemented
  void operator=(const Self &);  //purposely not implemented
  void TimeStepTestStability();
  double GeneralizedDiffCurve();
  double EdgeWeightedController(InputPixelType idxValue, InputPixelType centerValue);
  double meanNeighbors(NeighborIteratorType neighbors);
  InputPixelType up,down,left,right;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkAnisotropicAnomalousDiffusionImageFilter.hxx"
#endif

#endif
