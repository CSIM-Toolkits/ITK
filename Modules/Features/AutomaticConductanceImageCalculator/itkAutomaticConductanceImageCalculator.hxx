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
#ifndef itkAutomaticConductanceImageCalculator_hxx
#define itkAutomaticConductanceImageCalculator_hxx

#include "itkAutomaticConductanceImageCalculator.h"
#include "itkNumericTraits.h"

//Canny noise estimator
#include <itkGradientMagnitudeImageFilter.h>
#include <itkStatisticsImageFilter.h>
#include <itkImageToHistogramFilter.h>

//MAD
#include <itkMedianImageFilter.h>
#include <itkSubtractImageFilter.h>
#include <itkAbsImageFilter.h>

//Morphological
#include "itkGrayscaleMorphologicalOpeningImageFilter.h"
#include "itkGrayscaleMorphologicalClosingImageFilter.h"
#include "itkBinaryBallStructuringElement.h"

namespace itk
{
/**
 * Constructor
 */
template< typename TInputImage >
AutomaticConductanceImageCalculator< TInputImage >
::AutomaticConductanceImageCalculator()
{
    m_Image = TInputImage::New();
    m_Kappa = static_cast<PixelType>(1);
    m_OptimizationMethod = CANNY;
    m_RegionSetByUser = false;
}

/**
 * Compute the conductance (m_Kappa) of m_Image
 */
template< typename TInputImage >
void
AutomaticConductanceImageCalculator< TInputImage >
::Compute(void)
{
    if ( !m_RegionSetByUser )
    {
        m_Region = m_Image->GetRequestedRegion();
    }

    switch (m_OptimizationMethod) {
    case CANNY:
        ComputeCanny(m_Image);
        break;
    case MAD:
        ComputeMAD(m_Image);
        break;
    case MORPHOLOGICAL:
        ComputeMorphological(m_Image);
    default:
        break;
    }
}


template< typename TInputImage >
void
AutomaticConductanceImageCalculator< TInputImage >
::ComputeCanny(ImageConstPointer image)
{
    //Gradient magnitude
    typedef itk::GradientMagnitudeImageFilter<TInputImage, TInputImage>   GradientFilter;
    typename GradientFilter::Pointer gradient = GradientFilter::New();
    gradient->SetInput(image);

    //Gradient magnitude histogram
    const unsigned int MeasurementVectorSize = 1; // Grayscale
    const unsigned int binsPerDimension = 256;

    typedef itk::Statistics::ImageToHistogramFilter< TInputImage >        ImageToHistogramFilterType;
    typename ImageToHistogramFilterType::HistogramType::MeasurementVectorType     lowerBound(binsPerDimension);
    lowerBound.Fill(0);

    typename ImageToHistogramFilterType::HistogramType::MeasurementVectorType     upperBound(binsPerDimension);
    upperBound.Fill(255);

    typename ImageToHistogramFilterType::HistogramType::SizeType                  size(MeasurementVectorSize);
    size.Fill(binsPerDimension);

    typename ImageToHistogramFilterType::Pointer gradToHist =  ImageToHistogramFilterType::New();
    gradToHist->SetInput( gradient->GetOutput() );
    gradToHist->SetHistogramBinMinimum( lowerBound );
    gradToHist->SetHistogramBinMaximum( upperBound );
    gradToHist->SetHistogramSize( size );
    gradToHist->Update();

    typename ImageToHistogramFilterType::HistogramType* histogramG = gradToHist->GetOutput();

    double cumulativeSum = 0.0, total = 0.0;
    for (int idx = 0; idx < histogramG->GetSize()[0]; ++idx) {
        total+=histogramG->GetFrequency(idx);
    }

    for (int idx = 0; idx < histogramG->GetSize()[0]; ++idx) {
        cumulativeSum+=histogramG->GetFrequency(idx);
        double percentage = cumulativeSum/total;
        if (percentage >= 0.9) {
            m_Kappa = static_cast<PixelType>(histogramG->GetBinMinFromValue(0,histogramG->GetMeasurement(idx,0)));
            break;
        }
    }
}

template< typename TInputImage >
void
AutomaticConductanceImageCalculator< TInputImage >
::ComputeMAD(ImageConstPointer image)
{
    // K = 1.4826 median( |gradI - median(gradI)| )
    //Gradient magnitude
    typedef itk::GradientMagnitudeImageFilter<TInputImage, TInputImage>   GradientFilter;
    typename GradientFilter::Pointer gradient = GradientFilter::New();
    gradient->SetInput(image);


    //median(|gradI|)
    typedef itk::MedianImageFilter<TInputImage, TInputImage>    MedianFilter;
    typename MedianFilter::Pointer medianGrad = MedianFilter::New();
    medianGrad->SetInput(gradient->GetOutput());

    //Median deviation (|gradI - median(gradI)|)
    typedef itk::SubtractImageFilter<TInputImage, TInputImage>  SubtractFilter;
    typename SubtractFilter::Pointer deviation = SubtractFilter::New();
    deviation->SetInput1(gradient->GetOutput());
    deviation->SetInput2(medianGrad->GetOutput());


    typedef itk::AbsImageFilter<TInputImage, TInputImage>       AbsFilter;
    typename AbsFilter::Pointer abs = AbsFilter::New();
    abs->SetInput(deviation->GetOutput());

    typename MedianFilter::Pointer MAD = MedianFilter::New();
    MAD->SetInput(abs->GetOutput());

    typedef itk::StatisticsImageFilter<TInputImage>     StatisticsFilter;
    typename StatisticsFilter::Pointer statImage = StatisticsFilter::New();
    statImage->SetInput(MAD->GetOutput());
    statImage->Update();

    m_Kappa = 1.4826 * statImage->GetMean();
}

template< typename TInputImage >
void
AutomaticConductanceImageCalculator< TInputImage >
::ComputeMorphological(ImageConstPointer image)
{
    // K = average(opening(I,s)) - average(closing(I,s))
    // structuring element
    typedef itk::BinaryBallStructuringElement< PixelType, InputImageDimension  > StructuringElementType;
    // define the opening and closing types
    typedef itk::GrayscaleMorphologicalOpeningImageFilter<TInputImage, TInputImage, StructuringElementType >  OpeningFilterType;
    typedef itk::GrayscaleMorphologicalClosingImageFilter<TInputImage, TInputImage, StructuringElementType >  ClosingFilterType;

    // Create structuring element
    StructuringElementType  structuringElement;
    structuringElement.SetRadius( 1 );
    structuringElement.CreateStructuringElement();

    // Create the opening closing filters
    typename OpeningFilterType::Pointer  opening  = OpeningFilterType::New();
    typename ClosingFilterType::Pointer  closing  = ClosingFilterType::New();
    // Setup the opening and closing methods
    opening->SetKernel(  structuringElement );
    closing->SetKernel(  structuringElement );

    // creation of the pipeline. The enhancement operation is given by:
    // Original Image + Top Hat Image - Bottom Hat Image
    opening->SetInput( image );
    closing->SetInput( image );

    typedef itk::StatisticsImageFilter<TInputImage> StatisticsFilter;
    typename StatisticsFilter::Pointer statOpen = StatisticsFilter::New();
    typename StatisticsFilter::Pointer statClose = StatisticsFilter::New();
    statOpen->SetInput(opening->GetOutput());
    statClose->SetInput(closing->GetOutput());
    statOpen->Update();
    statClose->Update();

    m_Kappa= std::abs(statOpen->GetMean() - statClose->GetMean());
}

template< typename TInputImage >
void
AutomaticConductanceImageCalculator< TInputImage >
::SetRegion(const RegionType & region)
{
    m_Region = region;
    m_RegionSetByUser = true;
}

template< typename TInputImage >
void
AutomaticConductanceImageCalculator< TInputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);

      os << indent << "Kappa: "
         << static_cast< typename NumericTraits< PixelType >::PrintType >( m_Kappa )
         << std::endl;
    itkPrintSelfObjectMacro( Image );
    os << indent << "Region: " << std::endl;
    m_Region.Print( os, indent.GetNextIndent() );
    os << indent << "Region set by User: " << m_RegionSetByUser << std::endl;
}
} // end namespace itk

#endif
