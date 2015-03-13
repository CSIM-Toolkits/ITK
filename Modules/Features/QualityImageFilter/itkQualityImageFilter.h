#ifndef __itkQualityImageFilter_h
#define __itkQualityImageFilter_h
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkNumericTraits.h"

namespace itk
{

template< typename TImage >
class ITK_EXPORT QualityImageFilter:
        public ImageToImageFilter< TImage, TImage >
{
public:
    /** Extract dimension from inputs images, where it is assumed there are with the same type. */
    itkStaticConstMacro(TImageDimension, unsigned int,
                        TImage::ImageDimension);

    /** Standard class typedefs. */
    typedef QualityImageFilter                                 Self;
    typedef ImageToImageFilter< TImage, TImage >                  Superclass;
    typedef SmartPointer< Self >                                  Pointer;
    typedef SmartPointer< const Self >                            ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(QualityImageFilter, ImageToImageFilter)

    /** Set the Reference image. */
    void SetReferenceImage(const TImage *refImage);

    /** Set the Image for comparison. */
    void SetCompareImage(const TImage *compareImage);

    double SNR();
    double RMSE();
    double MAE();
    double SSIM();

protected:
    QualityImageFilter();
    virtual ~QualityImageFilter() {}
    void GenerateData();

private:
    QualityImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    double stdCorrelation();
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkQualityImageFilter.hxx"
#endif

#endif
