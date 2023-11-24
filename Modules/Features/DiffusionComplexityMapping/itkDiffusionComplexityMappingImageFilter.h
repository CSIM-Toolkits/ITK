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
#ifndef __itkDiffusionComplexityMappingImageFilter_h
#define __itkDiffusionComplexityMappingImageFilter_h
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

template< typename TInputImage, typename TOutputImage=Image<typename NumericTraits<typename TInputImage::ValueType>::ValueType, 3>, typename TInputMask=Image<unsigned char, 3> >
class ITK_EXPORT DiffusionComplexityMappingImageFilter:
        public ImageToImageFilter< TInputImage, TOutputImage >
{
public:
    /** Extract dimension from inputs images, where it is assumed there are with the same type. */
    itkStaticConstMacro(InputImageDimension, unsigned int,
                        TInputImage::ImageDimension);
    itkStaticConstMacro(OutputImageDimension, unsigned int,
                        TOutputImage::ImageDimension);

    /** Convenient typedefs for simplifying declarations. */
    typedef TInputImage     InputImageType;
    typedef TOutputImage    OutputImageType;
    typedef TInputMask      MaskImageType;

    /** Standard class typedefs. */
    typedef DiffusionComplexityMappingImageFilter          Self;
    typedef ImageToImageFilter< TInputImage, TOutputImage >       Superclass;
    typedef SmartPointer< Self >                                  Pointer;
    typedef SmartPointer< const Self >                            ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(DiffusionComplexityMappingImageFilter, ImageToImageFilter)

    typedef typename InputImageType::PixelType                  InputPixelType;
    typedef typename OutputImageType::PixelType                 OutputPixelType;
    typedef typename MaskImageType::PixelType                   MaskPixelType;
    typedef typename itk::Statistics::Histogram< OutputPixelType, itk::Statistics::DenseFrequencyContainer2 > HistogramType;

    void SetInputImage(const TInputImage* image);
    void SetDiffusionSpace(const TInputMask* mask);

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

    /** Set the disequilibrium function used in the complexity calculation. */
    itkSetMacro(DisequilibriumFunction, unsigned char)

    itkGetMacro(QValue, float)
    itkGetMacro(HistogramBins, unsigned int)
    itkGetMacro(UseManualNumberOfBins, bool)
    itkGetMacro(DisequilibriumFunction, unsigned char)
    itkGetMacro(DebugMode, bool)

#ifdef ITK_USE_CONCEPT_CHECKING
    // Begin concept checking
    itkConceptMacro( InputHasNumericTraitsCheck,
                     ( Concept::HasNumericTraits< InputPixelType > ) );
    itkConceptMacro( SameDimensionCheck,
                     ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );
#endif

    enum DisequilibriumFunction {
        EUCLIDEAN=1,
        WOOTER=2,
        KULLBACK_LEIBER=3,
        JENSEN=4
    };

protected:
    DiffusionComplexityMappingImageFilter();
    virtual ~DiffusionComplexityMappingImageFilter() {}
    virtual void GenerateOutputInformation(void) ITK_OVERRIDE;

    typename TInputImage::ConstPointer GetDWIImage();
    typename TInputMask::ConstPointer GetDiffusionSpace();
    void GenerateData();
private:
    DiffusionComplexityMappingImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    void createDiffusionSpace(typename InputImageType::Pointer diffImg, typename InputImageType::ConstPointer inputImg, std::vector<unsigned int> gradientsList);
    void getSpaceMaximumMinimumDiffusion(typename InputImageType::Pointer diffImg,typename MaskImageType::Pointer mask, OutputPixelType& maximum, OutputPixelType& minimum);
    void createDiffusionWeightedValues(typename InputImageType::Pointer diffAcquitions, typename InputImageType::Pointer diffImg, unsigned int numberOfGradients, unsigned int b0);
    void calculatesEntropyMapping(typename OutputImageType::Pointer output, typename InputImageType::Pointer diffImg, typename MaskImageType::Pointer mask, OutputPixelType max, OutputPixelType min);
    void calculatesDisequilibriumMapping(typename OutputImageType::Pointer output, typename InputImageType::Pointer diffImg, typename MaskImageType::Pointer mask, OutputPixelType max, OutputPixelType min, unsigned int numberGradients);
    unsigned int histogramBinByDiffusionVarianceRegression(unsigned int numberOfGradients, double max, double min);

    float m_QValue;
    unsigned int m_HistogramBins;
    bool m_UseManualNumberOfBins, m_DebugMode;
    unsigned char m_DisequilibriumFunction;
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkDiffusionComplexityMappingImageFilter.hxx"
#endif

#endif

