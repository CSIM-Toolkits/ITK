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
#ifndef __itkGlobalLogisticContrastEnhancementImageFilter_h
#define __itkGlobalLogisticContrastEnhancementImageFilter_h
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkNumericTraits.h"

namespace itk
{

template< typename TInputImage , typename TOutputImage>
class ITK_EXPORT GlobalLogisticContrastEnhancementImageFilter:
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
    typedef GlobalLogisticContrastEnhancementImageFilter          Self;
    typedef ImageToImageFilter< TInputImage, TOutputImage >       Superclass;
    typedef SmartPointer< Self >                                  Pointer;
    typedef SmartPointer< const Self >                            ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(GlobalLogisticContrastEnhancementImageFilter, ImageToImageFilter)

    typedef typename InputImageType::PixelType                 InputPixelType;
    typedef typename OutputImageType::PixelType                OutputPixelType;

    /** Set the maximum output. */
    itkSetMacro(MaximumOutput, double)

    /** Set the minimum output. */
    itkSetMacro(MinimumOutput, double)

    /** Set the object area. */
    itkSetMacro(FlipObjectArea, bool)
    itkBooleanMacro(FlipObjectArea)

    /** Set the higher cutting for logistic curve. */
    itkSetMacro(HigherCut, double)

    /** Set the lower cutting for logistic curve. */
    itkSetMacro(LowerCut, double)

    itkGetMacro(FlipObjectArea, bool)
    itkGetMacro(Alpha, double)
    itkGetMacro(Beta, double)
    itkGetMacro(MaximumOutput, double)
    itkGetMacro(MinimumOutput, double)

#ifdef ITK_USE_CONCEPT_CHECKING
    // Begin concept checking
    itkConceptMacro( InputHasNumericTraitsCheck,
                     ( Concept::HasNumericTraits< InputPixelType > ) );
    itkConceptMacro( SameDimensionCheck,
                     ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );
#endif

protected:
    GlobalLogisticContrastEnhancementImageFilter();
    virtual ~GlobalLogisticContrastEnhancementImageFilter() {}
    bool m_FlipObjectArea; // Use it if you want to invert the increasing of contrast of the image
    double m_Alpha, m_Beta; // Used as the logistic curve parameters
    double m_MaximumOutput, m_MinimumOutput; // Used to adjust the output values
    double m_HigherCut, m_LowerCut; // Used to remove the image outliers (default: the Brain Extraction Tool-FSL criteria)

    void GenerateData();
private:
    GlobalLogisticContrastEnhancementImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGlobalLogisticContrastEnhancementImageFilter.hxx"
#endif

#endif
