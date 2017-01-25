/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef itkDiffusionEdgeOptimizationImageCalculator_h
#define itkDiffusionEdgeOptimizationImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"

namespace itk
{
/** \class DiffusionEdgeOptimizationImageCalculator
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
 * \wikiexample{ImageProcessing/DiffusionEdgeOptimizationImageCalculator,Find the minimum and maximum value (and the position of the value) in an image}
 * \wikiexample{Developer/OilPaintingImageFilter,Multi-threaded oil painting image filter}
 * \endwiki
 */
template< typename TInputImage >
class DiffusionEdgeOptimizationImageCalculator:public Object
{
public:
    /** Extract dimension from input and output image. */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);

    /** Standard class typedefs. */
    typedef DiffusionEdgeOptimizationImageCalculator Self;
    typedef Object                        Superclass;
    typedef SmartPointer< Self >          Pointer;
    typedef SmartPointer< const Self >    ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(DiffusionEdgeOptimizationImageCalculator, Object);

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

    /** . */
    itkSetMacro(OptimizationMethod, unsigned char)

    /** Compute the minimum and maximum values of intensity of the input image. */
    void Compute();

    /** Return the minimum intensity value. */
    itkGetConstMacro(Kappa, PixelType);

    /** . */
    itkGetMacro(OptimizationMethod, unsigned char)

    /** Set the region over which the values will be computed */
    void SetRegion(const RegionType & region);

    enum OptimizationMethod {
        CANNY=1,
        MAD=2,
        MORPHOLOGICAL=3,
    };

protected:
    DiffusionEdgeOptimizationImageCalculator();
    virtual ~DiffusionEdgeOptimizationImageCalculator() {}
    virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:
    DiffusionEdgeOptimizationImageCalculator(const Self &) ITK_DELETE_FUNCTION;
    void operator=(const Self &) ITK_DELETE_FUNCTION;

    void ComputeCanny(ImageConstPointer image);
    void ComputeMAD(ImageConstPointer image);
    void ComputeMorphological(ImageConstPointer image);

    PixelType         m_Kappa;
    ImageConstPointer m_Image;
    unsigned char m_OptimizationMethod;

    RegionType m_Region;
    bool       m_RegionSetByUser;
};
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionEdgeOptimizationImageCalculator.hxx"
#endif

#endif /* itkDiffusionEdgeOptimizationImageCalculator_h */
