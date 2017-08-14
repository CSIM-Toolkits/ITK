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
#ifndef __itkDiffusionEntropyMappingImageFilter_hxx
#define __itkDiffusionEntropyMappingImageFilter_hxx
#include "itkDiffusionEntropyMappingImageFilter.h"

#include <itkNrrdImageIO.h>
#include <itkImageRegionIterator.h>
#include <itkDiffusionTensor3DReconstructionImageFilter.h>
#include <itkMetaDataObject.h>

#include <iostream>
using namespace std;

namespace itk
{
template< typename TInput, typename TOutput >
DiffusionEntropyMappingImageFilter< TInput, TOutput >
::DiffusionEntropyMappingImageFilter()
{
    m_QValue=1.0;
    m_HistogramBins=1;
}

template< typename TInput, typename TOutput >
void
DiffusionEntropyMappingImageFilter< TInput, TOutput >
::GenerateOutputInformation(void){
    OutputImageType *output = this->GetOutput();
    output->CopyInformation(this->GetInput());
    output->SetRegions(this->GetInput()->GetRequestedRegion());
    output->SetNumberOfComponentsPerPixel( 1 );

    output->Allocate();
    output->FillBuffer(0);
}

template< typename TInput, typename TOutput >
void
DiffusionEntropyMappingImageFilter< TInput, TOutput >
::GenerateData()
{
    unsigned int numberOfImages = 0;
    unsigned int numberOfGradientImages = 0;
    bool readb0 = false;
    double b0 = 0;
    std::vector<unsigned int> gradients;
    typename InputImageType::ConstPointer input = this->GetInput();
    typename OutputImageType::Pointer output = this->GetOutput();

    typedef itk::DiffusionTensor3DReconstructionImageFilter<OutputPixelType, OutputPixelType, double > TensorReconstructionImageFilterType;
    // -------------------------------------------------------------------------
    // Parse the Nrrd headers to get the B value and the gradient directions used
    // for diffusion weighting.
    //
    // The Nrrd headers should look like :
    // The tags specify the B value and the gradient directions. If gradient
    // directions are (0,0,0), it indicates that it is a reference image.
    //
    // DWMRI_b-value:=800
    // DWMRI_gradient_0000:= 0 0 0
    // DWMRI_gradient_0001:=-1.000000       0.000000        0.000000
    // DWMRI_gradient_0002:=-0.166000       0.986000        0.000000
    // DWMRI_gradient_0003:=0.110000        0.664000        0.740000
    // ...
    //
    itk::MetaDataDictionary imgMetaDictionary = input->GetMetaDataDictionary();
    std::vector<std::string> imgMetaKeys = imgMetaDictionary.GetKeys();
    std::vector<std::string>::const_iterator itKey = imgMetaKeys.begin();
    std::string metaString;
    typename TensorReconstructionImageFilterType::GradientDirectionType vect3d;
    typename TensorReconstructionImageFilterType::GradientDirectionContainerType::Pointer DiffusionVectors = TensorReconstructionImageFilterType::GradientDirectionContainerType::New();

    for (; itKey != imgMetaKeys.end(); ++itKey)
    {
        double x,y,z;
        itk::ExposeMetaData<std::string> (imgMetaDictionary, *itKey, metaString);
        if (itKey->find("DWMRI_gradient") != std::string::npos)
        {
            std::cout << *itKey << " ---> " << metaString << std::endl;
            sscanf(metaString.c_str(), "%lf %lf %lf\n", &x, &y, &z);
            vect3d[0] = x; vect3d[1] = y; vect3d[2] = z;
            DiffusionVectors->InsertElement( numberOfImages, vect3d );
            ++numberOfImages;
            // If the direction is 0.0, this is a reference image
            if (vect3d[0] == 0.0 &&
                    vect3d[1] == 0.0 &&
                    vect3d[2] == 0.0)
            {
                gradients.push_back(0);
                continue;
            }
            //            Add a new gradient volume
            gradients.push_back(1);
            ++numberOfGradientImages;
        }
        else if (itKey->find("DWMRI_b-value") != std::string::npos)
        {
            std::cout << *itKey << " ---> " << metaString << std::endl;
            readb0 = true;
            b0 = atof(metaString.c_str());
        }
    }
    std::cout << "Number of gradient images: "
              << numberOfGradientImages
              << " and Number of reference images: "
              << numberOfImages - numberOfGradientImages
              << std::endl;
    if(!readb0)
    {
        std::cerr << "BValue not specified in header file" << std::endl;
    }

    typedef itk::ImageRegionConstIterator<InputImageType>   RegionConstIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType>       RegionIteratorType;
    RegionConstIteratorType   diffusionIterator(input,input->GetRequestedRegion());
    RegionIteratorType        outputIterator(output, output->GetRequestedRegion());

    //    Iterate over the diffusion space to get the voxel data vector
    diffusionIterator.GoToBegin();
    outputIterator.GoToBegin();
    while (!diffusionIterator.IsAtEnd()) {
        //        Extract the vector data from voxel.
        //        This iterative process runs through each voxel of the input image
        //        where the diffusion gradients are allocated. With each voxel values
        //        trends, it is possible to reconstruct a local histogram and then calculate
        //        the local entropy.

        //Finding input minimum and maximum values
        OutputPixelType minimumInputValue = NumericTraits<OutputPixelType>::max();
        for (int min = 0; min < input->GetNumberOfComponentsPerPixel(); ++min) {
            if (diffusionIterator.Get()[min]<minimumInputValue){
                minimumInputValue=diffusionIterator.Get()[min];
            }
        }

        OutputPixelType maximumInputValue = NumericTraits<OutputPixelType>::min();
        for (int max = 0; max < input->GetNumberOfComponentsPerPixel(); ++max) {
            if (diffusionIterator.Get()[max]>maximumInputValue){
                maximumInputValue=diffusionIterator.Get()[max];
            }
        }

        typedef itk::Statistics::Histogram< typename OutputImageType::PixelType, itk::Statistics::DenseFrequencyContainer2 > HistogramType;
        typename HistogramType::Pointer histogram = HistogramType::New();
        typename HistogramType::SizeType size(1);

        //Setting the number of bins depending on a bin function.
        unsigned int binsPerDimension = automaticHistogramBinCalculation(input->GetNumberOfComponentsPerPixel());
        size.Fill(binsPerDimension);

        typename HistogramType::MeasurementVectorType lowerBound;
        lowerBound.SetSize(binsPerDimension);
        lowerBound.Fill(minimumInputValue);

        typename HistogramType::MeasurementVectorType upperBound;
        upperBound.SetSize(binsPerDimension);
        upperBound.Fill(maximumInputValue);
        histogram->SetMeasurementVectorSize(1);
        histogram->Initialize(size, lowerBound, upperBound );

        //Mounting normalized histogram
        //Step 1: Normalizing local data
        double Sum = 0;
        for (int i = 0; i < input->GetNumberOfComponentsPerPixel(); ++i) {
            if (gradients[i]!=0) {
                Sum += static_cast<double>(diffusionIterator.Get()[i]);
            }
        }
        typename HistogramType::IndexType index(1);
        typename HistogramType::MeasurementVectorType mv(1);

        //Step 2: Recalculating histogram
        for (int j = 0; j < input->GetNumberOfComponentsPerPixel(); ++j) {
            if (gradients[j]!=0) {
                mv[0]=diffusionIterator.Get()[j];
                histogram->GetIndex(mv,index);
                histogram->IncreaseFrequencyOfIndex(index, 1);
            }
        }

        //Calculate entropy
        OutputPixelType entropy = 0, p;
        for (int i = 0; i < histogram->GetSize()[0]; ++i) {
            if (m_QValue != 1.0) {
                p = static_cast<double>(histogram->GetFrequency(i))/static_cast<double>(histogram->GetTotalFrequency());
                if( p > NumericTraits<double>::min() )
                {
                    entropy += pow<double>(p,m_QValue);
                }
            }else{
                p = static_cast<double>(histogram->GetFrequency(i))/static_cast<double>(histogram->GetTotalFrequency());
                if( p > NumericTraits<double>::min() )
                {
                    entropy += - (p)
                            * (std::log( p ) / std::log( 2.0 ));
                }
            }
        }

        if (m_QValue != 1.0) {
            outputIterator.Set((1.0 - entropy)/(m_QValue - 1.0));
        }else{
            outputIterator.Set(entropy);
        }

        ++diffusionIterator;
        ++outputIterator;
    }

}

template< typename TInput, typename TOutput >
unsigned int
DiffusionEntropyMappingImageFilter< TInput, TOutput >
::automaticHistogramBinCalculation(unsigned int n)
{
    switch (m_HistogramBins) {
    case SQUAREROOT:
        return static_cast<unsigned int>(sqrt(n));
        break;
    case STURGES:
        return static_cast<unsigned int>(log2(n)+1);
        break;
    case RICE:
        return static_cast<unsigned int>(2*n^(1/3));
        break;
    default:
        break;
    }
}

} // end namespace itk

#endif
