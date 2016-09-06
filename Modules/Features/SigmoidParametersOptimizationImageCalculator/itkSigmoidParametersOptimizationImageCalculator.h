#ifndef __itkSigmoidParametersOptimizationImageCalculator_h
#define __itkSigmoidParametersOptimizationImageCalculator_h
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkNumericTraits.h"
#include "itkImageToHistogramFilter.h"

namespace itk
{

template< typename TImage >
class ITK_EXPORT SigmoidParametersOptimizationImageCalculator:
        public ImageToImageFilter< TImage, TImage >
{
public:
    /** Extract dimension from inputs images, where it is assumed there are with the same type. */
    itkStaticConstMacro(TImageDimension, unsigned int,
                        TImage::ImageDimension);

    /** Standard class typedefs. */
    typedef SigmoidParametersOptimizationImageCalculator          Self;
    typedef ImageToImageFilter< TImage, TImage >                  Superclass;
    typedef SmartPointer< Self >                                  Pointer;
    typedef SmartPointer< const Self >                            ConstPointer;

    typedef typename itk::Statistics::ImageToHistogramFilter< TImage >::HistogramPointer     ImageToHistogramType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(SigmoidParametersOptimizationImageCalculator, ImageToImageFilter)

    /** Set the Input image. */
    void SetInput(const TImage *input);

    /** Get the optimum sigmoid function parameters. */
    itkGetMacro(OptimumAlpha, double)
    itkGetMacro(OptimumBeta, double)

    /** Set the Energy kernel size. */
    itkSetMacro(KernelRadius, int)

    /** Set the Beta range. */
    itkSetMacro(MaximumBeta, double)
    itkSetMacro(MinimumBeta, double)

    /** Set the Alpha range. */
    itkSetMacro(MaximumAlpha, double)
    itkSetMacro(MinimumAlpha, double)

    /** Set the Histogram bins. */
    itkSetMacro(NumberOfBins, unsigned int)

    /** Set the number of trial to search the optimum parameters. */
    itkSetMacro(NumberOfTrials, unsigned int)

    itkGetMacro(KernelRadius, int)
    itkGetMacro(MaximumBeta, double)
    itkGetMacro(MinimumBeta, double)
    itkGetMacro(MaximumAlpha, double)
    itkGetMacro(MinimumAlpha, double)
    itkGetMacro(NumberOfBins, unsigned int)
    itkGetMacro(NumberOfTrials, unsigned int)

protected:
    SigmoidParametersOptimizationImageCalculator();
    virtual ~SigmoidParametersOptimizationImageCalculator() {}
    double m_OptimumBeta;
    double m_OptimumAlpha;
    int m_KernelRadius;
    double m_MaximumBeta;
    double m_MinimumBeta;
    double m_MaximumAlpha;
    double m_MinimumAlpha;
    unsigned int m_NumberOfBins;
    int m_NumberOfTrials;
    void GenerateData();

private:
    SigmoidParametersOptimizationImageCalculator(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    double sigmoid(double x, double alpha, double beta);
    double noiseFunction(ImageToHistogramType histogram, double alpha, double beta);
    double varianceFunction(ImageToHistogramType histogram, double alpha, double beta);
    double nonLinearFunction(ImageToHistogramType histogram, double alpha, double beta);
    //Implement Aj e Bj
    double Aj(ImageToHistogramType histogram,double alpha, double beta);
    double Bj(ImageToHistogramType histogram,double alpha, double beta);
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSigmoidParametersOptimizationImageCalculator.hxx"
#endif

#endif
