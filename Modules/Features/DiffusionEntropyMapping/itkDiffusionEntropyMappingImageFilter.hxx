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
template< typename TInput, typename TOutput, typename TMask >
DiffusionEntropyMappingImageFilter< TInput, TOutput, TMask >
::DiffusionEntropyMappingImageFilter()
{
    m_QValue=1.0;
    m_HistogramBins=2;
    m_UseManualNumberOfBins=false;
    this->SetNumberOfRequiredInputs(1);
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionEntropyMappingImageFilter<TInput, TOutput, TMask>::SetInputImage(const TInput* image)
{
    this->SetNthInput(0, const_cast<TInput*>(image));
}

template< typename TInput, typename TOutput, typename TMask>
void
DiffusionEntropyMappingImageFilter<TInput, TOutput, TMask>::SetDiffusionSpace(const TMask* mask)
{
    this->SetNthInput(1, const_cast<TMask*>(mask));
}

template< typename TInput, typename TOutput, typename TMask >
typename TInput::ConstPointer
DiffusionEntropyMappingImageFilter<TInput, TOutput, TMask>::GetDWIImage()
{
    return static_cast< const TInput * >
            ( this->ProcessObject::GetInput(0) );
}

template< typename TInput, typename TOutput, typename TMask >
typename TMask::ConstPointer
DiffusionEntropyMappingImageFilter<TInput, TOutput, TMask>::GetDiffusionSpace()
{
    return static_cast< const TMask * >
            ( this->ProcessObject::GetInput(1) );
}


template< typename TInput, typename TOutput, typename TMask >
void
DiffusionEntropyMappingImageFilter< TInput, TOutput, TMask >
::GenerateOutputInformation(void){
    OutputImageType *output = this->GetOutput();
    output->CopyInformation(this->GetDWIImage());
    output->SetRegions(this->GetDWIImage()->GetRequestedRegion());
    output->SetNumberOfComponentsPerPixel( 1 );

    output->Allocate();
    output->FillBuffer(0);
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionEntropyMappingImageFilter< TInput, TOutput, TMask >
::GenerateData()
{
    unsigned int numberOfImages = 0;
    unsigned int numberOfGradientImages = 0;
    bool readb0 = false;
    double b0 = 0;
    std::vector<unsigned int> gradients;
    typename InputImageType::ConstPointer input = this->GetDWIImage();
    typename OutputImageType::Pointer output = this->GetOutput();
    typename MaskImageType::Pointer usedMaskSpace;
    usedMaskSpace = MaskImageType::New();
    usedMaskSpace->CopyInformation(output);
    usedMaskSpace->SetRegions(output->GetRequestedRegion());
    usedMaskSpace->Allocate();

    if (this->GetNumberOfIndexedInputs()==1) {
        usedMaskSpace->FillBuffer(1);
    }else{
        typedef itk::ImageRegionIterator<MaskImageType>             InputRegionIteratorType;
        typedef itk::ImageRegionConstIterator<MaskImageType>        InputRegionConstIteratorType;
        InputRegionConstIteratorType   inputMaskIt(this->GetDiffusionSpace(), this->GetDiffusionSpace()->GetRequestedRegion());
        InputRegionIteratorType        usedMaskIt(usedMaskSpace, usedMaskSpace->GetRequestedRegion());

        inputMaskIt.GoToBegin();
        usedMaskIt.GoToBegin();
        while (!inputMaskIt.IsAtEnd()) {
            usedMaskIt.Set(inputMaskIt.Get());
            ++inputMaskIt;
            ++usedMaskIt;
        }
    }

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
            if (m_DebugMode) {
                std::cout << *itKey << " ---> " << metaString << std::endl;
            }
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
            if (m_DebugMode) {
                std::cout << *itKey << " ---> " << metaString << std::endl;
            }
            readb0 = true;
            b0 = atof(metaString.c_str());
        }
    }
    if (m_DebugMode) {
        std::cout << "Number of gradient images: "
                  << numberOfGradientImages
                  << " and Number of reference images: "
                  << numberOfImages - numberOfGradientImages
                  << std::endl;
    }
    if(!readb0)
    {
        std::cerr << "BValue not specified in header file" << std::endl;
    }

    //First: take the mean S0 image if there are more than one non-diffusion volume in the image
    //Creating the new input image with only one b0 volume (the mean of the N non-diffusion images)
    typename InputImageType::Pointer diffusionAcquisitionImage = InputImageType::New();
    diffusionAcquisitionImage->CopyInformation(input);
    diffusionAcquisitionImage->SetRegions(input->GetRequestedRegion());
    diffusionAcquisitionImage->SetVectorLength( numberOfGradientImages + 1 );
    diffusionAcquisitionImage->Allocate();

    createDiffusionSpace(diffusionAcquisitionImage, input, gradients);

    //Secondly: Calculates the diffusion weights for each direction (D = (ln(Si) - ln(S0))/b )
    typename InputImageType::Pointer diffusionImage = InputImageType::New();
    diffusionImage->CopyInformation(input);
    diffusionImage->SetRegions(input->GetRequestedRegion());
    diffusionImage->SetVectorLength( numberOfGradientImages );
    diffusionImage->Allocate();

    createDiffusionWeightedValues(diffusionAcquisitionImage, diffusionImage, numberOfGradientImages, b0);

    //Finding input minimum and maximum values only on the gradient volumes and also construct the global probability distribution
    OutputPixelType minimumInputValue = NumericTraits<OutputPixelType>::max();
    OutputPixelType maximumInputValue = NumericTraits<OutputPixelType>::min();
    getSpaceMaximumMinimumDiffusion(diffusionImage, usedMaskSpace, maximumInputValue, minimumInputValue);

    if (m_DebugMode) {
        std::cout<<"Diffusion space range (max: "<<maximumInputValue<<" - min: "<<minimumInputValue<<")"<<std::endl;
    }

    //Setting the number of bins depending on a bin function.
    if (m_UseManualNumberOfBins) {
        if (m_DebugMode) {
            cout<<"Manual number of bins: "<<m_HistogramBins<<endl;
        }
    }else{
        m_HistogramBins = numberOfGradientImages;
        if (m_DebugMode) {
            cout<<"Automatic number of bins: "<<m_HistogramBins<<endl;
        }
    }

    calculatesEntropyMapping(output, diffusionImage, usedMaskSpace, maximumInputValue, minimumInputValue);
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionEntropyMappingImageFilter< TInput, TOutput, TMask >
::createDiffusionSpace(typename InputImageType::Pointer diffImg, typename InputImageType::ConstPointer inputImg, std::vector<unsigned int> gradientsList)
{
    unsigned int numberOfGradientImages=0, numberOfReferenceImages=0;
    typedef itk::ImageRegionConstIterator<InputImageType>   InputRegionConstIteratorType;
    typedef itk::ImageRegionIterator<InputImageType>        InputRegionIteratorType;
    InputRegionConstIteratorType   inputIterator(inputImg,inputImg->GetRequestedRegion());
    InputRegionIteratorType    rawDiffIt(diffImg, diffImg->GetRequestedRegion());

    //Detect the number of gradients images
    for (int i = 0; i < gradientsList.size(); ++i) {
        if (gradientsList[i]!=0) {
            numberOfGradientImages++;
        }
    }

    numberOfReferenceImages=inputImg->GetNumberOfComponentsPerPixel() - numberOfGradientImages;

    typedef itk::VariableLengthVector<OutputPixelType> DiffusionVectorType;

    rawDiffIt.GoToBegin();
    inputIterator.GoToBegin();
    while (!inputIterator.IsAtEnd()) {
        DiffusionVectorType       weights;
        weights.SetSize(numberOfGradientImages + 1);

        //Capturing the mean non-diffusion values
        OutputPixelType meanNonDiff=static_cast<OutputPixelType>(0);
        for (unsigned int n = 0; n < inputImg->GetNumberOfComponentsPerPixel(); ++n) {
            if (gradientsList[n]==0){
                meanNonDiff+=inputIterator.Get()[n];
            }
        }
        meanNonDiff/=numberOfReferenceImages;
        weights[0]=meanNonDiff;

        unsigned int index=1;
        for (unsigned int n = 0; n < inputImg->GetNumberOfComponentsPerPixel(); ++n) {
            if (gradientsList[n]!=0){
                weights[index]=inputIterator.Get()[n];
                ++index;
            }
        }
        rawDiffIt.Set(weights);

        ++rawDiffIt;
        ++inputIterator;
    }
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionEntropyMappingImageFilter< TInput, TOutput, TMask >
::getSpaceMaximumMinimumDiffusion(typename InputImageType::Pointer diffImg, typename MaskImageType::Pointer mask, OutputPixelType& maximum, OutputPixelType& minimum)
{
    typedef itk::ImageRegionIterator<InputImageType>        InputRegionIteratorType;
    typedef itk::ImageRegionConstIterator<MaskImageType>        InputMaskRegionConstIteratorType;
    InputRegionIteratorType    diffIt(diffImg, diffImg->GetRequestedRegion());
    InputMaskRegionConstIteratorType    maskIt(mask, mask->GetRequestedRegion());

    diffIt.GoToBegin();
    maskIt.GoToBegin();
    while (!diffIt.IsAtEnd()) {
        if (maskIt.Get()!=static_cast<MaskPixelType>(0)) {
            for (int min = 0; min < diffImg->GetNumberOfComponentsPerPixel(); ++min) {
                if (diffIt.Get()[min]<minimum && diffIt.Get()[min]>0){
                    minimum=diffIt.Get()[min];
                }
            }

            for (int max = 0; max < diffImg->GetNumberOfComponentsPerPixel(); ++max) {
                if (diffIt.Get()[max]>maximum && !isinf(diffIt.Get()[max])){
                    maximum=diffIt.Get()[max];
                }
            }
        }
        ++diffIt;
        ++maskIt;
    }
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionEntropyMappingImageFilter< TInput, TOutput, TMask >
::createDiffusionWeightedValues(typename InputImageType::Pointer diffAcquitions, typename InputImageType::Pointer diffImg, unsigned int numberOfGradientImages, unsigned int b0)
{
    typedef itk::ImageRegionIterator<InputImageType>        InputRegionIteratorType;
    InputRegionIteratorType    rawDiffIt(diffAcquitions, diffAcquitions->GetRequestedRegion());
    InputRegionIteratorType    diffIt(diffImg, diffImg->GetRequestedRegion());

    typedef itk::VariableLengthVector<OutputPixelType> DiffusionVectorType;

    diffIt.GoToBegin();
    rawDiffIt.GoToBegin();
    while (!rawDiffIt.IsAtEnd()) {
        DiffusionVectorType       diffValues;
        diffValues.SetSize(numberOfGradientImages);

        unsigned int index=0;
        for (unsigned int n = 1; n < diffAcquitions->GetNumberOfComponentsPerPixel(); ++n) {
            if (rawDiffIt.Get()[0]!=0 || rawDiffIt.Get()[n]!=0) {
                diffValues[index]=std::abs((std::log(rawDiffIt.Get()[n]) - std::log(rawDiffIt.Get()[0]))/static_cast<OutputPixelType>(b0));
            }else{
                diffValues[index]=0;
            }
            ++index;
        }
        diffIt.Set(diffValues);

        ++diffIt;
        ++rawDiffIt;
    }


    //    typedef itk::ImageFileWriter<InputImageType>    WriterType;
    //    typename WriterType::Pointer w = WriterType::New();
    //    w->SetInput(diffImg);
    //    w->SetFileName("/home/antonio/Downloads/DiffusionEntropyMappingImageFilter/dwi_f_hard_data.nii.gz");
    //    w->Update();
}

template< typename TInput, typename TOutput, typename TMask >
void
DiffusionEntropyMappingImageFilter< TInput, TOutput, TMask >
::calculatesEntropyMapping(typename OutputImageType::Pointer output, typename InputImageType::Pointer diffImg, typename MaskImageType::Pointer mask, OutputPixelType max, OutputPixelType min)
{
    //    Iterate over the entire diffusion space to get the information contained in each voxel
    typedef itk::ImageRegionIterator<InputImageType>        InputRegionIteratorType;
    typedef itk::ImageRegionIterator<OutputImageType>       OutputRegionIteratorType;
    typedef itk::ImageRegionIterator<MaskImageType>         MaskRegionIteratorType;
    InputRegionIteratorType     diffIt(diffImg, diffImg->GetRequestedRegion());
    OutputRegionIteratorType    outputIterator(output, output->GetRequestedRegion());
    MaskRegionIteratorType      maskIterator(mask, mask->GetRequestedRegion());
    diffIt.GoToBegin();
    outputIterator.GoToBegin();
    maskIterator.GoToBegin();

    //Normalize the entropy value regarding the equiprobable states ln(N)
    double maxEntropy = std::log(m_HistogramBins);
    while (!outputIterator.IsAtEnd()) {
        if (maskIterator.Get()>static_cast<MaskPixelType>(0)) {
            //        Extract the vector data from voxel.
            //        This iterative process runs through each voxel of the input image
            //        where the diffusion gradients are allocated. With each voxel values
            //        trends, it is possible to reconstruct a local histogram and then calculate
            //        the local entropy.
            //Create the PDF for i-eth voxel
            typename HistogramType::Pointer diffusionProbability = HistogramType::New();
            typename HistogramType::SizeType size(1);

            size.Fill(m_HistogramBins);

            typename HistogramType::MeasurementVectorType lowerBound;
            lowerBound.SetSize(m_HistogramBins);
            lowerBound.Fill(min);

            typename HistogramType::MeasurementVectorType upperBound;
            upperBound.SetSize(m_HistogramBins);
            upperBound.Fill(max);
            diffusionProbability->SetMeasurementVectorSize(1);
            diffusionProbability->Initialize(size, lowerBound, upperBound );

            typename HistogramType::IndexType index(1);
            typename HistogramType::MeasurementVectorType mv(1);

            for (int j = 0; j < diffImg->GetNumberOfComponentsPerPixel(); ++j) {
                if ( diffIt.Get()[j]>static_cast<OutputPixelType>(0) && !isinf(diffIt.Get()[j]) ) {
                    mv[0]=diffIt.Get()[j];
                    diffusionProbability->GetIndex(mv,index);
                    diffusionProbability->IncreaseFrequencyOfIndex(index, 1);
                }
            }

            //Calculate entropy
            OutputPixelType entropy = 0.0, p;
            for (int j = 0; j < diffImg->GetNumberOfComponentsPerPixel(); ++j) {
                p = static_cast<double>(diffusionProbability->GetFrequency(j))/static_cast<double>(diffusionProbability->GetTotalFrequency());
                if (p==0.0) {
                    continue;
                }else{
                    if (m_QValue != 1.0) {
                        if( p > NumericTraits<double>::min() )
                        {
                            entropy += pow<double>(p,m_QValue);
                        }
                    }else{
                        if( p > NumericTraits<double>::min() )
                        {
                            entropy += - (p)*(std::log( p ));
                        }
                    }
                }
            }

            if (m_QValue != 1.0) {
                outputIterator.Set(((1.0 - entropy)/(m_QValue - 1.0))/maxEntropy);
            }else{
                outputIterator.Set(entropy/maxEntropy);
            }
        }else{
            outputIterator.Set(0);
        }

        ++diffIt;
        ++outputIterator;
        ++maskIterator;
    }
}

} // end namespace itk

#endif
