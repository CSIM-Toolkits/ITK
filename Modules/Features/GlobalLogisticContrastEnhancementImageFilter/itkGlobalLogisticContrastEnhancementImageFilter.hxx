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
#ifndef __itkGlobalLogisticContrastEnhancementImageFilter_hxx
#define __itkGlobalLogisticContrastEnhancementImageFilter_hxx
#include "itkGlobalLogisticContrastEnhancementImageFilter.h"

#include <itkImageToHistogramFilter.h>
#include <itkImageRegionConstIterator.h>
#include <itkImageRegionIterator.h>
#include <itkStatisticsImageFilter.h>
#include <itkSigmoidImageFilter.h>
#include <itkMultiplyImageFilter.h>

#include <math.h>

using namespace std;

namespace itk
{
template< typename TInput, typename TOutput>
GlobalLogisticContrastEnhancementImageFilter< TInput, TOutput >
::GlobalLogisticContrastEnhancementImageFilter()
{
    this->m_FlipObjectArea=false;
    this->m_Alpha=0.0;
    this->m_Beta=0.0;
    this->m_MaximumOutput=1.0;
    this->m_MinimumOutput=0.0;
    this->m_LowerCut=0.02; // based on the Brain Extraction Tool (FSL) criteria
    this->m_HigherCut=0.98;
}

template< typename TInput, typename TOutput >
void
GlobalLogisticContrastEnhancementImageFilter< TInput, TOutput >
::GenerateData()
{
    //Input image
    typename InputImageType::ConstPointer input = this->GetInput();
    //Output image
    typename OutputImageType::Pointer output = this->GetOutput();
    output->SetRegions(input->GetBufferedRegion());
    output->Allocate();

    //Image statistics
    typedef itk::StatisticsImageFilter<InputImageType> StatisticsImageFilterType;
    typename StatisticsImageFilterType::Pointer inputStatistics = StatisticsImageFilterType::New ();
    inputStatistics->SetInput(input);
    inputStatistics->Update();

    typedef itk::Statistics::ImageToHistogramFilter< InputImageType >   HistogramFilterType;
    typename HistogramFilterType::Pointer histogramFilter = HistogramFilterType::New();

    typedef typename HistogramFilterType::HistogramSizeType   SizeType;
    SizeType size( 1 );
    size[0] = std::sqrt(inputStatistics->GetMaximum() - inputStatistics->GetMinimum());

    histogramFilter->SetHistogramSize( size );
    histogramFilter->SetMarginalScale( 10.0 );

    typename HistogramFilterType::HistogramMeasurementVectorType lowerBound( 1 );
    typename HistogramFilterType::HistogramMeasurementVectorType upperBound( 1 );
    lowerBound[0] = ((inputStatistics->GetMaximum()-inputStatistics->GetMinimum())/size[0]);
    upperBound[0] = inputStatistics->GetMaximum();
    histogramFilter->SetHistogramBinMinimum( lowerBound );
    histogramFilter->SetHistogramBinMaximum( upperBound );

    histogramFilter->SetInput(  input  );
    histogramFilter->Update();

    typedef typename HistogramFilterType::HistogramType  HistogramType;
    const HistogramType *histogram = histogramFilter->GetOutput();

    //Setting the logistic curve
    double max_logistic = histogram->Quantile(0,m_HigherCut), min_logistic = histogram->Quantile(0,m_LowerCut);
    m_Beta = (max_logistic - min_logistic)/2.0; // centered the logistic curve in the maximum range of the input image
    m_Alpha = inputStatistics->GetSigma(); // Assuming the image standard deviation in order to get the maximum of information in the input image

    if (m_FlipObjectArea) {
        m_Alpha=(-1)*m_Alpha;
    }

    //Apply sigmoid on input image
    typedef itk::SigmoidImageFilter<InputImageType, OutputImageType> SigmoidFilterType;
    typename SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
    sigmoid->SetInput(input);
    sigmoid->SetOutputMinimum(m_MinimumOutput);
    sigmoid->SetOutputMaximum(m_MaximumOutput);
    sigmoid->SetAlpha(m_Alpha);
    sigmoid->SetBeta(m_Beta);
    sigmoid->Update();

    //Applying the contrast enhancement function on the input image
    typedef itk::MultiplyImageFilter<InputImageType, OutputImageType>  MultiplyFilterType;
    typename MultiplyFilterType::Pointer mul = MultiplyFilterType::New();
    mul->SetInput1(input);
    mul->SetInput2(sigmoid->GetOutput());
    mul->Update();

    itk::ImageRegionConstIterator<TInput> contrastEnhancementIterator(mul->GetOutput(), mul->GetOutput()->GetBufferedRegion());
    itk::ImageRegionIterator<TOutput> outputIterator(output, output->GetBufferedRegion());

    contrastEnhancementIterator.IsAtBegin();
    outputIterator.IsAtBegin();
    while (!contrastEnhancementIterator.IsAtEnd()) {
        outputIterator.Set(contrastEnhancementIterator.Get());
        ++contrastEnhancementIterator;
        ++outputIterator;
    }

}

} // end namespace itk

#endif
