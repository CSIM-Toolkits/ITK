
#ifndef __itkGeneralizedEntropyThresholdImageCalculator_h
#define __itkGeneralizedEntropyThresholdImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkNumericTraits.h"

namespace itk
{

/** \class GeneralizedEntropyThresholdImageCalculator
 * \brief Computes the MaxEntropy's threshold for an image.
 * 
 *
 * Implements Kapur-Sahoo-Wong (Maximum Entropy) thresholding method
 * Kapur J.N., Sahoo P.K., and Wong A.K.C. (1985) "A New Method for
 * Gray-Level Picture Thresholding Using the Entropy of the Histogram"
 * Graphical Models and Image Processing, 29(3): 273-285
 * M. Emre Celebi
 * 06.15.2007
 * Ported to ImageJ plugin by G.Landini from E Celebi's fourier_0.8 routines
 *
 *
 * Ported from the ImageJ implementation. 
 *
 * This class is templated over the input image type.
 * \author Richard Beare
 * \warning This method assumes that the input image consists of scalar pixel
 * types.
 *
 * \ingroup Operators
 */
template <class TInputImage>
class ITK_EXPORT GeneralizedEntropyThresholdImageCalculator : public Object
{
public:
  /** Standard class typedefs. */
  typedef GeneralizedEntropyThresholdImageCalculator Self;
  typedef Object                       Superclass;
  typedef SmartPointer<Self>           Pointer;
  typedef SmartPointer<const Self>     ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro(GeneralizedEntropyThresholdImageCalculator, Object);

  /** Type definition for the input image. */
  typedef TInputImage  ImageType;

  /** Pointer type for the image. */
  typedef typename TInputImage::Pointer  ImagePointer;
  
  /** Const Pointer type for the image. */
  typedef typename TInputImage::ConstPointer ImageConstPointer;

  /** Type definition for the input image pixel type. */
  typedef typename TInputImage::PixelType PixelType;
  
  /** Type definition for the input image region type. */
  typedef typename TInputImage::RegionType RegionType;
  
  /** Set the input image. */
  itkSetConstObjectMacro(Image,ImageType);

  /** Compute the MaxEntropy's threshold for the input image. */
  void Compute(void);

  /** Return the MaxEntropy's threshold value. */
  itkGetConstMacro(Threshold,PixelType);

  /** Set the q parameter. */
  itkSetMacro(Q,double);
  
  /** Set/Get the number of histogram bins. Default is 128. */
  itkSetClampMacro( NumberOfHistogramBins, unsigned long, 1, 
                    NumericTraits<unsigned long>::max() );
  itkGetConstMacro( NumberOfHistogramBins, unsigned long );


  /** Set the region over which the values will be computed */
  void SetRegion( const RegionType & region );

protected:
  GeneralizedEntropyThresholdImageCalculator();
  virtual ~GeneralizedEntropyThresholdImageCalculator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

private:
  GeneralizedEntropyThresholdImageCalculator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
  
  PixelType            m_Threshold;
  double               m_Q;
  unsigned long        m_NumberOfHistogramBins;
  ImageConstPointer    m_Image;
  RegionType           m_Region;
  bool                 m_RegionSetByUser;

};

} // end namespace itk


#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGeneralizedEntropyThresholdImageCalculator.hxx"
#endif

#endif
