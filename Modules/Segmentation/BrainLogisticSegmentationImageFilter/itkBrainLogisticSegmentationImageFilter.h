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
#ifndef __itkBrainLogisticSegmentationImageFilter_h
#define __itkBrainLogisticSegmentationImageFilter_h
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
class ITK_EXPORT BrainLogisticSegmentationImageFilter:
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
    typedef BrainLogisticSegmentationImageFilter          Self;
    typedef ImageToImageFilter< TInputImage, TOutputImage >       Superclass;
    typedef SmartPointer< Self >                                  Pointer;
    typedef SmartPointer< const Self >                            ConstPointer;

    /** Method for creation through the object factory. */
    itkNewMacro(Self)

    /** Run-time type information (and related methods). */
    itkTypeMacro(BrainLogisticSegmentationImageFilter, ImageToImageFilter)

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
    BrainLogisticSegmentationImageFilter();
    virtual ~BrainLogisticSegmentationImageFilter() {}
    bool m_ManualTolerance, m_DebugMode, m_UseManualNumberOfBins;
    char m_Tolerance;
    unsigned int m_NumberOfBins, m_NumberOfTissues; //TODO Fazer o GetSet do m_NumberOfTissues

    virtual void GenerateOutputInformation(void) ITK_OVERRIDE;
    void GenerateData();
private:
    BrainLogisticSegmentationImageFilter(const Self &); //purposely not implemented
    void operator=(const Self &);  //purposely not implemented
    void checkTolerance(char tolerance);
};

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkBrainLogisticSegmentationImageFilter.hxx"
#endif

#endif
