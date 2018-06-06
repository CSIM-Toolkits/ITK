#ifndef __itkSampEn2DImageCalculator_h
#define __itkSampEn2DImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"
#include "itkImageRegionConstIterator.h"
#include "itkStatisticsImageFilter.h"
#include "itkNumericTraits.h"

#include <stdlib.h>
#include <math.h>
#include <ctime>

namespace itk
{
/** \class SampEn2DImageCalculator
 *  \brief Computes the minimum and the maximum intensity values of
 *         an image.
 *
 * This calculator computes the minimum and the maximum intensity values of
 * an image.  It is templated over input image type.  If only Maximum or
 * Minimum value is needed, just call ComputeMaximum() (ComputeMinimum())
 * otherwise Compute() will compute both.
 *
 * \ingroup Operators
 * \ingroup ITKCommon
 *
 * \wiki
 * \wikiexample{ImageProcessing/MinimumMaximumImageCalculator,Find the minimum and maximum value (and the position of the value) in an image}
 * \wikiexample{Developer/OilPaintingImageFilter,Multi-threaded oil painting image filter}
 * \endwiki
 */
template< typename TInputImage >
class SampEn2DImageCalculator:public Object
{
public:
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);

    /** Standard class typedefs. */
    typedef SampEn2DImageCalculator Self;
    typedef Object                        Superclass;
    typedef SmartPointer< Self >          Pointer;
    typedef SmartPointer< const Self >    ConstPointer;

    /** Type definition for the output image entropy value. This is the type which the image will be cast. */
    typedef bool RFlagType;
    typedef double DoublePixelType;
    typedef unsigned int MParameterValueType;
    typedef double RParameterValueType;
    typedef unsigned int DParameterValueType;
    typedef double BParameterValueType;
    typedef ImageRegionConstIterator<TInputImage>       ConstRegionIteratorType;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(SampEn2DImageCalculator, Object);

    /** Type definition for the input image. */
    typedef TInputImage ImageType;

    /** Pointer type for the image. */
    typedef typename TInputImage::Pointer ImagePointer;

    /** Const Pointer type for the image. */
    typedef typename TInputImage::ConstPointer ImageConstPointer;

    /** Type definition for the input image pixel type. */
    typedef typename TInputImage::PixelType PixelType;


    /** Type definition for the input image index type. */
    typedef typename TInputImage::IndexType IndexType;

    /** Type definition for the input image region type. */
    typedef typename TInputImage::RegionType RegionType;

    /** Set the input image. */
    itkSetConstObjectMacro(Image, ImageType);

    /** Set the RParameterAsPercentage parameter value. */
    itkSetMacro(UseRParameterAsPercentage, RFlagType)
    itkBooleanMacro(UseRParameterAsPercentage)

    /** Set the M parameter value. */
    itkSetMacro(M, MParameterValueType);

    /** Set the R parameter value. */
    itkSetMacro(R, RParameterValueType);

    /** Set the D (delay) parameter value. */
    itkSetMacro(D, DParameterValueType);

    /** Set the B (background value) parameter value. */
    itkSetMacro(BGV, BParameterValueType);

    /** Get the M parameter value. */
    itkGetMacro(M, MParameterValueType);

    /** Get the R parameter value. */
    itkGetMacro(R, RParameterValueType);

    /** Get the D parameter value. */
    itkGetMacro(D, DParameterValueType);

    /** Get the B parameter value. */
    itkGetMacro(BGV, BParameterValueType);

    /** Compute the two-dimensional sample entropy value of the input image. */
    void ComputeEntropy(void);

    /** Return the entropy value. */
    itkGetConstMacro(Entropy, DoublePixelType);

    /** Set the region over which the values will be computed */
    void SetRegion(const RegionType & region);


protected:
    SampEn2DImageCalculator();
    virtual ~SampEn2DImageCalculator() {}
    virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;


private:
    SampEn2DImageCalculator(const Self &); //purposely not implemented
    void operator=(const Self &);                //purposely not implemented

    RFlagType           m_UseRParameterAsPercentage;
    DoublePixelType     m_Entropy;
    MParameterValueType m_M;
    RParameterValueType m_R;
    DParameterValueType m_D;
    BParameterValueType m_BGV;
    ImageConstPointer   m_Image;

    /** Test to check the M, R, D and B parameters values */
    void ParametersCertification(); 

    RegionType m_Region;
    bool       m_RegionSetByUser;
    size_t        m_Nx,m_Ny;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSampEn2DImageCalculator.hxx"
#endif

#endif /* __itkSampEn2DImageCalculator_h */
