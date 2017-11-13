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

    /** Set threshold method. */
    itkSetMacro(ThresholdMethod, unsigned char)

    itkGetMacro(FlipObjectArea, bool)
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
        ISODATA=6,
        INTERMODES=7,
    };

protected:
    LogisticContrastEnhancementImageFilter();
    virtual ~LogisticContrastEnhancementImageFilter() {}
    bool m_FlipObjectArea;
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
