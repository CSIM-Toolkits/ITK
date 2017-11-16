/*
   Copyright 2017 Antonio Carlos da Silva Senra Filho and Fabricio Henrique Simozo

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
#ifndef __itkDiffusionEntropyMappingImageFilter_h
#define __itkDiffusionEntropyMappingImageFilter_h
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkNumericTraits.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "itkHistogram.h"
#include "itkDenseFrequencyContainer2.h"
#include <math.h>

namespace itk
{

template< typename TInputImage, typename TOutputImage=Image<typename NumericTraits<typename TInputImage::ValueType>::ValueType, 3> >
class ITK_EXPORT DiffusionEntropyMappingImageFilter:
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
    typedef DiffusionEntropyMappingImageFilter          Self;
    typedef ImageToImageFilter< TInputImage, TOutputImage >       Superclass;
    typedef SmartPointer< Self >                                  Pointer;
    typedef SmartPointer< const Self >                            ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(DiffusionEntropyMappingImageFilter, ImageToImageFilter)

    typedef typename InputImageType::PixelType                 InputPixelType;
    typedef typename OutputImageType::PixelType                OutputPixelType;

    /** Set the q value used in the entropy calculation. */
    itkSetMacro(QValue, float)

    /** Choose if define the number of bins manually. */
    itkBooleanMacro(UseManualNumberOfBins)
    itkSetMacro(UseManualNumberOfBins, bool)

    /** Set the number of bins used in the entropy calculation. */
    itkSetMacro(HistogramBins, unsigned int)

    /** Debug mode is used to inform some messages in the standard output. */
    itkBooleanMacro(DebugMode)
    itkSetMacro(DebugMode, bool)

    itkGetMacro(QValue, float)
    itkGetMacro(HistogramBins, unsigned int)
    itkGetMacro(UseManualNumberOfBins, bool)
    itkGetMacro(DebugMode, bool)

#ifdef ITK_USE_CONCEPT_CHECKING
    // Begin concept checking
    itkConceptMacro( InputHasNumericTraitsCheck,
                     ( Concept::HasNumericTraits< InputPixelType > ) );
    itkConceptMacro( SameDimensionCheck,
                     ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );
#endif

    enum BinsEstimate {
        SQUAREROOT=1,
        STURGES=2,
        RICE=3
    };

protected:
    DiffusionEntropyMappingImageFilter();
    virtual ~DiffusionEntropyMappingImageFilter() {}
    virtual void GenerateOutputInformation(void) ITK_OVERRIDE;
    void GenerateData();
private:
    DiffusionEntropyMappingImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    unsigned int automaticHistogramBinCalculation(unsigned int n);
    float m_QValue;
    unsigned int m_HistogramBins;
    bool m_UseManualNumberOfBins, m_DebugMode;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionEntropyMappingImageFilter.hxx"
#endif

#endif
