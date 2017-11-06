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
#ifndef __itkMultipleLogisticClassificationImageFilter_hxx
#define __itkMultipleLogisticClassificationImageFilter_hxx
#include "itkMultipleLogisticClassificationImageFilter.h"




using namespace std;

namespace itk
{
template< typename TInput, typename TOutput>
MultipleLogisticClassificationImageFilter< TInput, TOutput >
::MultipleLogisticClassificationImageFilter()
{
    this->m_Tolerance=1;
    this->m_ManualTolerance=false;
    this->m_DebugMode=false;
    //    this->m_NumberOfBins=32;
    this->m_NumberOfTissues=3;
}

template< typename TInput, typename TOutput >
void
MultipleLogisticClassificationImageFilter< TInput, TOutput >
::GenerateOutputInformation(void){
    this->Superclass::GenerateOutputInformation();

    OutputImageType *output = this->GetOutput();
    output->SetNumberOfComponentsPerPixel( m_NumberOfTissues );
}

template< typename TInput, typename TOutput >
void
MultipleLogisticClassificationImageFilter< TInput, TOutput >
::GenerateData()
{
    checkTolerance(m_Tolerance);

    //Input image
    typename InputImageType::ConstPointer input = this->GetInput();
    //Output image
    typename OutputImageType::Pointer output = static_cast< OutputImageType * >(this->GetOutput());
    output->SetRegions(input->GetBufferedRegion());
    output->Allocate();

    //Mounting input histogram
    //Step 1: Getting input data
    typedef itk::ImageRegionConstIterator<InputImageType> ConstRegionIterator;
    ConstRegionIterator inputIt(input, input->GetRequestedRegion());
    inputIt.GoToBegin();

    //Step 2: Removing zeros
    inputIt.GoToBegin();
    std::vector<InputPixelType> data;
    while (!inputIt.IsAtEnd()) {
        if (inputIt.Get() > static_cast<InputPixelType>(0)) {
            data.push_back(inputIt.Get());
        }
        ++inputIt;
    }

    //Step 3:Finding minimum and maximum of the input image
    InputPixelType min=NumericTraits<InputPixelType>::max(),max=NumericTraits<InputPixelType>::min();
    for (int i = 0; i < data.size(); ++i) {
        if (data[i]<min) {
            min=data[i];
        }
        if (data[i]>max) {
            max=data[i];
        }
    }
    if (m_DebugMode) {
        cout<<"Input (min,max) values considered in the histogram: ("<<min<<","<<max<<")"<<endl;
    }

    //Step 4:Defining the number of bins from input values range
    std::vector<InputPixelType> peaks, valleys;
    //Setting the first values to zero
    for (int i = 0; i < m_NumberOfTissues+1; ++i) {
        peaks.push_back(0);
    }
    for (int i = 0; i < m_NumberOfTissues+1; ++i) {
        valleys.push_back(0);
    }


    bool doBinCorrection=true;
    int binCorrectionWeight = 1;
    while(doBinCorrection){
        if (!m_UseManualNumberOfBins) {
            m_NumberOfBins=(int)(sqrt(max-min))*binCorrectionWeight; //Using the rule of squares.
        }

        //Image histogram
        typedef itk::Statistics::Histogram< typename InputImageType::PixelType, itk::Statistics::DenseFrequencyContainer2 > HistogramType;
        typename HistogramType::Pointer histogram = HistogramType::New();
        typename HistogramType::SizeType size(1);

        unsigned int binsPerDimension = m_NumberOfBins;
        size.Fill(binsPerDimension);

        typename HistogramType::MeasurementVectorType lowerBound;
        lowerBound.SetSize(binsPerDimension);
        lowerBound.Fill(min);

        typename HistogramType::MeasurementVectorType upperBound;
        upperBound.SetSize(binsPerDimension);
        upperBound.Fill(max);
        histogram->SetMeasurementVectorSize(1);
        histogram->Initialize(size, lowerBound, upperBound );


        typename HistogramType::IndexType index(1);
        typename HistogramType::MeasurementVectorType mv(1);

        //Step 5: Recalculating histogram
        for (int j = 0; j < data.size(); ++j) {
            mv[0]=data[j];
            histogram->GetIndex(mv,index);
            histogram->IncreaseFrequencyOfIndex(index, 1);
        }


        //Step 6: Removing outliers in the input data. Procesure similar to what is done in FSL-BET brain extraction.
        vector<int> freq;
        for (int p = 0; p < histogram->Size(); ++p) {
            freq.push_back(histogram->GetFrequency(p));
        }

        vector<double> cdf;
        InputPixelType sum=0,cdf_sum=0;
        for (int p = 0; p < freq.size(); ++p) {
            sum+=freq[p];
        }

        for (int i = 0; i < freq.size(); ++i) {
            cdf_sum+=(double)freq[i]/(double)sum;
            cdf.push_back(cdf_sum);
        }

        for (int i = 0; i < cdf.size(); ++i) {
            //Correcting data set to remove CDF lower than 2% and higher than 98%.
            if (cdf[i]<0.02 || cdf[i]>0.99) {
                index(0)=i;
                histogram->SetFrequencyOfIndex(index, 0);
            }
        }

        //Step 7: Finding n peaks and (n-1) valleys in the corrected histogram



        //Calculating 1st derivate of the histogram
        std::vector<double> derivate;
        for (int p = 0; p < m_NumberOfBins - 1; ++p) {
            derivate.push_back(((double)histogram->GetFrequency(p+1)-(double)histogram->GetFrequency(p)));
        }
        //Getting the derivative zeros. The peaks are found from high to low voxel values.
        std::vector<unsigned int> max_index, min_index;
        for (int i = derivate.size(); i > 0; --i) {
            if (derivate[i]*derivate[i-1]<0 && derivate[i]<0) {
                max_index.push_back(i-1);
            }
            if (derivate[i]*derivate[i-1]<0 && derivate[i]>0) {
                min_index.push_back(i-1);
            }
        }

        //Getting the peaks from the histogram
        for (int n = 0; n < max_index.size(); ++n) {
            peaks[n]=histogram->GetBinMax(0,max_index[n]);
            m_HistogramPeaks.push_back(static_cast<double>(peaks[n]));
        }
        for (int n = 0; n < min_index.size(); ++n) {
            valleys[n]=histogram->GetBinMax(0,min_index[n]);
            m_HistogramValleys.push_back(static_cast<double>(valleys[n]));
        }

        if (m_DebugMode && binCorrectionWeight>1) {
            cout<<"Correcting number of bins: "<<binCorrectionWeight<<" times..."<<endl;
        }
        if (peaks[m_NumberOfTissues-1]!=0 || binCorrectionWeight>3) { //A stop criteria is needed in order to not increase the nuber of bins too much.
            doBinCorrection=false;
        }
        binCorrectionWeight++;
    }

    if (m_DebugMode) {
        cout<<"Number of bins: "<<m_NumberOfBins<<endl;
        cout<<"Histogram peaks: ";
        for (int n = 0; n < m_NumberOfTissues; ++n) {
            cout<<peaks[n]<<" ";
        }
        cout<<endl;
        cout<<"Histogram valleys: ";
        for (int n = 0; n < m_NumberOfTissues; ++n) {
            cout<<valleys[n]<<" ";
        }
        cout<<endl;
    }


    //Step 8:Set Beta
    /* The shiftting factor is set based on the peaks and valleys interval.
     * In this way, the central value of the logistic function will be set
     * on the local minimum value found in the image histogram.
     */
    InputPixelType beta=0.0, alpha=0.0;
    std::vector<typename InputImageType::Pointer> tissues;
    for (int n = 0; n < m_NumberOfTissues; ++n) {
        typename InputImageType::Pointer outputTissue = InputImageType::New();
        outputTissue->SetRegions(input->GetBufferedRegion());
        outputTissue->Allocate();
        tissues.push_back(outputTissue);
    }

    typedef itk::SigmoidImageFilter<InputImageType, InputImageType>                SigmoidFilterType;
    typedef itk::MultiplyImageFilter<InputImageType, InputImageType>                MultiplyFilterType;
    typename SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
    typename MultiplyFilterType::Pointer mul = MultiplyFilterType::New();

    typename InputImageType::Pointer i1 = InputImageType::New();
    typename InputImageType::Pointer i2 = InputImageType::New();
    i1->SetRegions(input->GetBufferedRegion());
    i1->Allocate();
    i2->SetRegions(input->GetBufferedRegion());
    i2->Allocate();

    /* Here it is where all the segmentation process is done.
     */
    for (int n = 0; n < m_NumberOfTissues; ++n) {
        /* The first and last peaks only need one logistic classification
         * , since they are on the extrema of the image
         */
        if (n==0) {
            beta = valleys[n];
            m_Betas.push_back(beta);
            //Adjusting automatic tolerance
            if (m_ManualTolerance) {
                alpha=((-1)*(peaks[n+1])+beta)/(log((100.0-static_cast<double>(m_Tolerance))/static_cast<double>(m_Tolerance)));
                m_Alphas.push_back(static_cast<double>(alpha));
            }else{
                alpha=((-1)*(peaks[n+1])+beta)/(log(0.99/(1.0-0.99)));
                m_Alphas.push_back(static_cast<double>(alpha));
            }

            //Apply logistic curve in the n-th tissue
            sigmoid->SetInput(input);
            sigmoid->SetOutputMinimum(0);
            sigmoid->SetOutputMaximum(1);
            sigmoid->SetAlpha(alpha);
            sigmoid->SetBeta(beta);
            sigmoid->Update();

            itk::ImageRegionConstIterator<TInput> sigmoidIt(sigmoid->GetOutput(), sigmoid->GetOutput()->GetBufferedRegion());
            itk::ImageRegionIterator<TInput> tissueIt(tissues[n], tissues[n]->GetBufferedRegion());

            //Calculation the low boundary image using the first setting of logistic classification
            sigmoidIt.IsAtBegin();
            tissueIt.IsAtBegin();
            while (!sigmoidIt.IsAtEnd()) {
                tissueIt.Set(sigmoidIt.Get());
                ++sigmoidIt;
                ++tissueIt;
            }
        }else if (n==m_NumberOfTissues - 1) {
            beta = valleys[n-1];
            m_Betas.push_back(beta);
            //Adjusting automatic tolerance
            if (m_ManualTolerance) {
                alpha=((-1)*(peaks[n-1])+beta)/(log((100.0-static_cast<double>(m_Tolerance))/static_cast<double>(m_Tolerance)));
                m_Alphas.push_back(static_cast<double>(alpha));
            }else{
                alpha=((-1)*(peaks[n-1])+beta)/(log(0.99/(1.0-0.99)));
                m_Alphas.push_back(static_cast<double>(alpha));
            }

            //Apply logistic curve in the n-th tissue
            sigmoid->SetInput(input);
            sigmoid->SetOutputMinimum(0);
            sigmoid->SetOutputMaximum(1);
            sigmoid->SetAlpha(alpha);
            sigmoid->SetBeta(beta);
            sigmoid->Update();

            itk::ImageRegionConstIterator<TInput> sigmoidIt(sigmoid->GetOutput(), sigmoid->GetOutput()->GetBufferedRegion());
            itk::ImageRegionConstIterator<TInput> inputIt(input, input->GetBufferedRegion());
            itk::ImageRegionIterator<TInput> tissueIt(tissues[n], tissues[n]->GetBufferedRegion());

            //Calculation the low boundary image using the first setting of logistic classification
            sigmoidIt.IsAtBegin();
            inputIt.IsAtBegin();
            tissueIt.IsAtBegin();
            while (!sigmoidIt.IsAtEnd()) {
                if (inputIt.Get()!=0) {
                    tissueIt.Set(sigmoidIt.Get());
                }
                ++sigmoidIt;
                ++inputIt;
                ++tissueIt;
            }
        }else{
            beta = valleys[n];
            m_Betas.push_back(beta);
            //Adjusting automatic tolerance
            if (m_ManualTolerance) {
                alpha=((-1)*(peaks[n+1])+beta)/(log((100.0-static_cast<double>(m_Tolerance))/static_cast<double>(m_Tolerance)));
                m_Alphas.push_back(static_cast<double>(alpha));
            }else{
                alpha=((-1)*(peaks[n+1])+beta)/(log(0.99/(1.0-0.99)));
                m_Alphas.push_back(static_cast<double>(alpha));
            }

            //Apply logistic curve in the n-th tissue
            sigmoid->SetInput(input);
            sigmoid->SetOutputMinimum(0);
            sigmoid->SetOutputMaximum(1);
            sigmoid->SetAlpha(alpha);
            sigmoid->SetBeta(beta);
            sigmoid->Update();

            itk::ImageRegionConstIterator<TInput> sigmoidIt(sigmoid->GetOutput(), sigmoid->GetOutput()->GetBufferedRegion());
            itk::ImageRegionIterator<TInput> i1It(i1, i1->GetBufferedRegion());

            //Calculation the low boundary image using the first setting of logistic classification
            sigmoidIt.IsAtBegin();
            i1It.IsAtBegin();
            while (!sigmoidIt.IsAtEnd()) {
                i1It.Set(sigmoidIt.Get());
                ++sigmoidIt;
                ++i1It;
            }

            beta=valleys[n-1];
            m_Betas.push_back(beta);
            if (m_ManualTolerance) {
                alpha=((-1)*(peaks[n-1])+beta)/(log((100.0-static_cast<double>(m_Tolerance))/static_cast<double>(m_Tolerance)));
                m_Alphas.push_back(static_cast<double>(alpha));
            }else{
                alpha=((-1)*(peaks[n-1])+beta)/(log(0.99/(1.0-0.99)));
                m_Alphas.push_back(static_cast<double>(alpha));
            }

            //Apply logistic curve in the n-th tissue
            sigmoid->SetInput(input);
            sigmoid->SetOutputMinimum(0);
            sigmoid->SetOutputMaximum(1);
            sigmoid->SetAlpha(alpha);
            sigmoid->SetBeta(beta);
            sigmoid->Update();

            itk::ImageRegionIterator<TInput> i2It(i2, i2->GetBufferedRegion());

            //Calculation the low boundary image using the first setting of logistic classification
            sigmoidIt.GoToBegin();
            i2It.GoToBegin();
            while (!sigmoidIt.IsAtEnd()) {
                i2It.Set(sigmoidIt.Get());
                ++sigmoidIt;
                ++i2It;
            }

            //Calculating fuzzy segmentation S = I1.I2
            mul->SetInput1(i1);
            mul->SetInput2(i2);
            mul->Update();

            itk::ImageRegionIterator<TInput> sIt(mul->GetOutput(), mul->GetOutput()->GetBufferedRegion());
            itk::ImageRegionIterator<TInput> tissueIt(tissues[n], tissues[n]->GetBufferedRegion());

            tissueIt.IsAtBegin();
            sIt.IsAtBegin();
            while (!sIt.IsAtEnd()) {
                tissueIt.Set(sIt.Get());
                ++tissueIt;
                ++sIt;
            }
        }
    }

    //Setting the output volume
    typedef itk::ComposeImageFilter<InputImageType>     ComposeVectorType;
    typename ComposeVectorType::Pointer outputVector = ComposeVectorType::New();
    for (int n = 0; n < m_NumberOfTissues; ++n) {
        outputVector->SetInput(n, tissues[n]);
    }
    outputVector->Update();

    itk::ImageRegionIterator<OutputImageType> outputIt(output, output->GetRequestedRegion());
    itk::ImageRegionIterator<OutputImageType> vectorIt(outputVector->GetOutput(), outputVector->GetOutput()->GetRequestedRegion());

    vectorIt.IsAtBegin();
    outputIt.IsAtBegin();
    while (!outputIt.IsAtEnd()) {
        outputIt.Set(vectorIt.Get());
        ++vectorIt;
        ++outputIt;
    }

}

template<typename TInput, typename TOutput>
void
MultipleLogisticClassificationImageFilter<TInput, TOutput>
::checkTolerance(char tolerance) {
    if (tolerance >= 50) {
        std::cout<<"Tolerance is out of bound. Set in percentual [0 to 50["<<std::endl;
        exit(EXIT_FAILURE);
    }
}

} // end namespace itk

#endif
