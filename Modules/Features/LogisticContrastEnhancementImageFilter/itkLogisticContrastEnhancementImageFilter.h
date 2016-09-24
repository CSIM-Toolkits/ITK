#ifndef __itkLogisticContrastEnhancementImageFilter_h
#define __itkLogisticContrastEnhancementImageFilter_h
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkNumericTraits.h"

namespace itk
{

template< typename TInputImage , typename TOutputImage>
class ITK_EXPORT LogisticContrastEnhancementImageFilter:
        public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
    /** Extract dimension from inputs images, where it is assumed there are with the same type. */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                        TOutputImage::ImageDimension);

    /** Convenient typedefs for simplifying declarations. */
    typedef TInputImage  InputImageType;
    typedef TOutputImage OutputImageType;


    /** Standard class typedefs. */
    typedef LogisticContrastEnhancementImageFilter          Self;
    typedef ImageToImageFilter< TInputImage, TOutputImage >       Superclass;
    typedef SmartPointer< Self >                                  Pointer;
    typedef SmartPointer< const Self >                            ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(LogisticContrastEnhancementImageFilter, ImageToImageFilter)

    typedef typename InputImageType::PixelType                 InputPixelType;
    typedef typename OutputImageType::PixelType                OutputPixelType;

    /** Set the maximum output. */
    itkSetMacro(MaximumOutput, double)

    /** Set the minimum output. */
    itkSetMacro(MinimumOutput, double)

    /** Set the Histogram bins. */
    itkSetMacro(NumberOfBins, unsigned int)

    /** Set the Tolerance. */
    itkSetMacro(Tolerance, char)

    /** Set the object area. */
    itkSetMacro(FlipObjectArea, bool)
    itkBooleanMacro(FlipObjectArea)

    /** Set if use manual tolerance. */
    itkSetMacro(ManualTolerance, bool)
    itkBooleanMacro(ManualTolerance)

    /** Set threshold method. */
    itkSetMacro(ThresholdMethod, unsigned char)

    itkGetMacro(FlipObjectArea, bool)
    itkGetMacro(ManualTolerance, bool)
    itkGetMacro(Alpha, double)
    itkGetMacro(Beta, double)
    itkGetMacro(MaximumOutput, double)
    itkGetMacro(MinimumOutput, double)
    itkGetMacro(NumberOfBins, unsigned int)
    itkGetMacro(Tolerance, char)
    itkGetMacro(ThresholdMethod, unsigned char)

#ifdef ITK_USE_CONCEPT_CHECKING
    // Begin concept checking
    itkConceptMacro( InputHasNumericTraitsCheck,
                     ( Concept::HasNumericTraits< InputPixelType > ) );
    itkConceptMacro( SameDimensionCheck,
                     ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );
#endif

    enum ThresholdMethod {
        MAXENTROPY=1,
        OTSU=2,
        RENYI=3,
        MOMENTS=4,
        YEN=5,
        ISODATA=6
    };

protected:
    LogisticContrastEnhancementImageFilter();
    virtual ~LogisticContrastEnhancementImageFilter() {}
    bool m_FlipObjectArea;
    bool m_ManualTolerance;
    char m_Tolerance;
    double m_Alpha;
    double m_Beta;
    double m_MaximumOutput;
    double m_MinimumOutput;
    unsigned int m_NumberOfBins;
    unsigned char m_ThresholdMethod;

    void GenerateData();
private:
//    const unsigned short MAX_TOLERANCE = 0.99;
    LogisticContrastEnhancementImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    double sigmoid(double x, double alpha, double beta);
    void checkTolerance(char tolerance);
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLogisticContrastEnhancementImageFilter.hxx"
#endif

#endif
