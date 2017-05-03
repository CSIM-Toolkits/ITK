/*
   Copyright 2016 Antonio Carlos da Silva Senra Filho

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
 */
#ifndef itkAutomaticConductanceImageCalculator_h
#define itkAutomaticConductanceImageCalculator_h

#include "itkObject.h"
#include "itkObjectFactory.h"

namespace itk
{
/** \class AutomaticConductanceImageCalculator
 *  \brief Computes the conductance parameter based on flux functions
 * for anisotropic diffusion filtering process.
 *
 * This calculator computes the conductance parameter, usually denoted as Kappa,
 * which uses a gradient flux function already presented in the following papers:
 *
 * 1. Perona, P. & Malik, J., 1990. Scale-space and edge detection using anisotropic diffusion.
 * IEEE Transactions on Pattern Analysis and Machine Intelligence, 12(7), pp.629–639.
 *
 * 2. Black, M.J. et al., 1998. Robust Anisotropic Diffusion.
 * IEEE Trans. Image Process., 7(3), p.421.
 *
 * 3. Voci, F. et al., 2004. Estimating the gradient threshold in the perona-malik equation.
 * IEEE Signal Processing Magazine, 23(3), pp.39–46.
 *
 * The input image is analysed for noise intensity and edge consistency, where a
 * determined conductance is estimated.
 */
template< typename TInputImage >
class AutomaticConductanceImageCalculator:public Object
{
public:
    /** Extract dimension from input and output image. */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);

    /** Standard class typedefs. */
    typedef AutomaticConductanceImageCalculator Self;
    typedef Object                        Superclass;
    typedef SmartPointer< Self >          Pointer;
    typedef SmartPointer< const Self >    ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self);

    /** Run-time type information (and related methods). */
    itkTypeMacro(AutomaticConductanceImageCalculator, Object);

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

    /** Compute the conductance based on the input image. */
    void Compute();

    /** Return the conductance intensity value. */
    itkGetConstMacro(Kappa, PixelType);

    /** Set the type of optimization function used. */
    itkGetMacro(OptimizationMethod, unsigned char)

    /** Set the region over which the values will be computed */
    void SetRegion(const RegionType & region);

    enum OptimizationMethod {
        CANNY=1,
        MAD=2,
        MORPHOLOGICAL=3,
    };

protected:
    AutomaticConductanceImageCalculator();
    virtual ~AutomaticConductanceImageCalculator() {}
    virtual void PrintSelf(std::ostream & os, Indent indent) const ITK_OVERRIDE;

private:
    AutomaticConductanceImageCalculator(const Self &) ITK_DELETE_FUNCTION;
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
#include "itkAutomaticConductanceImageCalculator.hxx"
#endif

#endif /* itkAutomaticConductanceImageCalculator_h */
