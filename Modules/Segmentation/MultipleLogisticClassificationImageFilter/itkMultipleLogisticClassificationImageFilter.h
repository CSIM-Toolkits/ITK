/*
   Copyright 2017 Antonio Carlos da Silva Senra Filho

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
#ifndef __itkMultipleLogisticClassificationImageFilter_h
#define __itkMultipleLogisticClassificationImageFilter_h
#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkNumericTraits.h"
#include <itkVariableLengthVector.h>
#include <itkComposeImageFilter.h>
#include <itkVectorImage.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkStatisticsImageFilter.h>
#include <itkHistogram.h>
#include <itkSigmoidImageFilter.h>
#include <itkMultiplyImageFilter.h>
#include <math.h>

namespace itk
{

template< typename TInputImage , typename TOutputImage=VectorImage<typename NumericTraits<typename TInputImage::ValueType>::ValueType, TInputImage::ImageDimension> >
class ITK_EXPORT MultipleLogisticClassificationImageFilter:
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
    typedef MultipleLogisticClassificationImageFilter          Self;
    typedef ImageToImageFilter< TInputImage, TOutputImage >       Superclass;
    typedef SmartPointer< Self >                                  Pointer;
    typedef SmartPointer< const Self >                            ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(MultipleLogisticClassificationImageFilter, ImageToImageFilter)

    typedef typename InputImageType::PixelType                 InputPixelType;
    typedef typename OutputImageType::PixelType                OutputPixelType;

    /** Set the Histogram bins. */
    itkSetMacro(NumberOfBins, unsigned int)

    /** Set the Number of Tissues to be segmented. */
    itkSetMacro(NumberOfTissues, unsigned int)

    /** Set the Tolerance. */
    itkSetMacro(Tolerance, char)

    /** Set if use manual tolerance. */
    itkSetMacro(ManualTolerance, bool)
    itkBooleanMacro(ManualTolerance)

    /** Set if use manual number of bins. */
    itkSetMacro(UseManualNumberOfBins, bool)
    itkBooleanMacro(UseManualNumberOfBins)

    /** Set if use debug mode. */
    itkSetMacro(DebugMode, bool)
    itkBooleanMacro(DebugMode)


    /** Get if use manual tolerance. */
    itkGetMacro(ManualTolerance, bool)

    /** Get if use manual number of bins. */
    itkGetMacro(UseManualNumberOfBins, bool)

    /** Get if use debug mode. */
    itkGetMacro(DebugMode, bool)

    /** Set the Histogram bins. */
    itkGetMacro(NumberOfBins, unsigned int)

    /** Set the Number of Tissues to be segmented. */
    itkGetMacro(NumberOfTissues, unsigned int)

    /** Get the Tolerance. */
    itkGetMacro(Tolerance, char)

    /** Get the vector of peaks found in the histogram analysis. */
    itkGetMacro(HistogramPeaks, std::vector<double>)

    /** Get the vector of valleys found in the histogram analysis. */
    itkGetMacro(HistogramValleys, std::vector<double>)

    /** Get the vector of alpha values estimated in each part of the histogram.
    The counter follows the same orientation that is given in the HistogramPeaks and HistogramValleys,
    i.e. from the higher to the lower values of the histogram.*/
    itkGetMacro(Alphas, std::vector<double>)

    /** Get the vector of beta values estimated in each part of the histogram.
    The counter follows the same orientation that is given in the HistogramPeaks and HistogramValleys,
    i.e. from the higher to the lower values of the histogram.*/
    itkGetMacro(Betas, std::vector<double>)

#ifdef ITK_USE_CONCEPT_CHECKING
    // Begin concept checking
    itkConceptMacro( InputHasNumericTraitsCheck,
                     ( Concept::HasNumericTraits< InputPixelType > ) );
    itkConceptMacro( SameDimensionCheck,
                     ( Concept::SameDimension< InputImageDimension, OutputImageDimension > ) );
    itkConceptMacro( InputCovertibleToOutputCheck,
                       ( Concept::Convertible< InputPixelType, typename NumericTraits<OutputPixelType>::ValueType > ) );
#endif

protected:
    MultipleLogisticClassificationImageFilter();
    virtual ~MultipleLogisticClassificationImageFilter() {}
    bool m_ManualTolerance, m_DebugMode, m_UseManualNumberOfBins;
    char m_Tolerance;
    std::vector<double> m_HistogramPeaks, m_HistogramValleys, m_Alphas, m_Betas;
    unsigned int m_NumberOfBins, m_NumberOfTissues; //TODO Fazer o GetSet do m_NumberOfTissues

    virtual void GenerateOutputInformation(void) ITK_OVERRIDE;
    void GenerateData();
private:
    MultipleLogisticClassificationImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    void checkTolerance(char tolerance);
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMultipleLogisticClassificationImageFilter.hxx"
#endif

#endif
