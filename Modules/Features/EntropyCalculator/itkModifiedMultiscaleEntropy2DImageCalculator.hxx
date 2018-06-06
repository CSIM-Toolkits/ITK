/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkModifiedMultiscaleEntropy2DImageCalculator_hxx
#define __itkModifiedMultiscaleEntropy2DImageCalculator_hxx

#include "itkModifiedMultiscaleEntropy2DImageCalculator.h"
#include "itkSampEn2DImageCalculator.h"
#include "itkImageRegionIterator.h"
#include "itkConstNeighborhoodIterator.h"

#include "itkImageFileWriter.h"

using namespace std;

namespace itk
{
/**
 * Constructor
 */
template< typename TInputImage >
ModifiedMultiscaleEntropy2DImageCalculator< TInputImage >
::ModifiedMultiscaleEntropy2DImageCalculator()
{
    m_Image = TInputImage::New();
    //    m_Entropy.clear();
    m_M = 1;
    m_R = 0.10;
    m_S = 10;
    m_BGV = NAN;
    m_UseRParameterAsPercentage=true;
    m_RegionSetByUser = false;
}

/**
 * Compute Entropy of m_Image
 */
template< typename TInputImage >
void
ModifiedMultiscaleEntropy2DImageCalculator< TInputImage >
::ComputeEntropy(void)
{
    if ( !m_RegionSetByUser )
    {
        m_Region = m_Image->GetRequestedRegion();
    }

    ParametersCertification();

    typename ImageType::SizeType radius;
    radius.Fill(static_cast<int>(m_M));

    typename ImageType::SizeType size = m_Region.GetSize();

    m_Nx = size[0]; //ncols
    m_Ny = size[1]; //nrows

    double image_matrix[m_Nx*m_Ny];


    //Assigning tolerance factor
    double tolerance;

    if(m_UseRParameterAsPercentage) {
        double sigma = 0.0, mean = 0.0;
        ConstRegionIteratorType    copyIt(m_Image, m_Region);
        copyIt.GoToBegin(); // mean
        int count=0;
        double N=0;

        //Calculating stantard deviation
        if(isnan(m_BGV)) { //whole image
            while (!copyIt.IsAtEnd()) {
                image_matrix[count++]=copyIt.Get();
                N++;
                mean+=copyIt.Get();
                ++copyIt;
            }
            mean/=N;

            copyIt.GoToBegin(); // std
            while (!copyIt.IsAtEnd()) {
                sigma+=pow(copyIt.Get()-mean,2.0);
                ++copyIt;
            }
        } else {
            //Only in the foreground area
            while (!copyIt.IsAtEnd()) {
                image_matrix[count++]=copyIt.Get();

                if (copyIt.Get()>m_BGV) {
                    N++;
                    mean+=copyIt.Get();
                }
                ++copyIt;
            }
            mean/=N;

            copyIt.GoToBegin(); // std
            while (!copyIt.IsAtEnd()) {
                if (copyIt.Get()>m_BGV) {
                    sigma+=pow(copyIt.Get()-mean,2.0);
                }
                ++copyIt;
            }
        }

        sigma/=(N-1);
        sigma=sqrt(sigma);
        //End standard deviation

        tolerance = m_R*sigma;
    } else {
        tolerance = m_R;
    }

    //Call SampeEn2D for S=1
    typedef itk::SampEn2DImageCalculator<TInputImage> SampEn2DType;
    typename SampEn2DType::Pointer sampEn2D = SampEn2DType::New();
    sampEn2D->UseRParameterAsPercentageOff();
    sampEn2D->SetImage(m_Image);
    sampEn2D->SetM(m_M);
    sampEn2D->SetR(tolerance);
    sampEn2D->SetD(1);
    sampEn2D->SetBGV(m_BGV);
    sampEn2D->ComputeEntropy();

    m_Entropy.push_back(sampEn2D->GetEntropy());

    //Loop for remaining S
    typedef itk::Image<double, InputImageDimension> ScaledImageType;
    typedef itk::SampEn2DImageCalculator<ScaledImageType> ScaledSampEn2DType;
    typename ScaledSampEn2DType::Pointer scaledSampEn2D = ScaledSampEn2DType::New();
    scaledSampEn2D->UseRParameterAsPercentageOff();
    scaledSampEn2D->SetM(m_M);
    scaledSampEn2D->SetR(tolerance);
    scaledSampEn2D->SetBGV(m_BGV);

    int nPixels, s, i, j, si, sj, w, h;
    double sum;
    for(s=2; s<=m_S; s++) {
        //Coarse-grain - ITK
        typename ScaledImageType::Pointer scaledImage = ScaledImageType::New();
        typename ImageType::SizeType scaledSize;
        scaledSize[0] = m_Nx - s + 1; //ncols
        scaledSize[1] = m_Ny - s + 1; //nrows

        typename TInputImage::RegionType scaledRegion;
        scaledRegion.SetSize(scaledSize);
        scaledRegion.SetIndex(m_Region.GetIndex());

        scaledImage->SetRegions(scaledRegion);
        scaledImage->Allocate();

        // Coarse-grain
        w = m_Nx-s+1;
        h = m_Ny-s+1;
        nPixels = s*s;

        double cg_image[w*h], checkBGV;
        bool calculateCG=true;
        for (i=0; i < h; i++) {
            for (j=0; j < w; j++) {
                sum = 0.0;

                for (si=0; si < s; si++) {
                    for (sj=0; sj < s; sj++) {
                        checkBGV=image_matrix[(i + si)*m_Nx+(j + sj)];
                        if (checkBGV>m_BGV) {
                            sum += checkBGV;
                        }else{
                            cg_image[i*w+j] = m_BGV;
                            sj=s;
                            si=s;
                            calculateCG=false;
                            break;
                        }
                    } // for si
                } // for sj

                cg_image[i*w+j] = (calculateCG)?sum/nPixels:m_BGV;
                calculateCG=true;
            } // for j
        } // for i

        typedef itk::ImageRegionIterator<ScaledImageType>   ImageRegionIteratorType;
        ImageRegionIteratorType scaledIt(scaledImage, scaledImage->GetRequestedRegion());
        scaledIt.GoToBegin();
        int count=0;
        while (!scaledIt.IsAtEnd()) {
            scaledIt.Set(cg_image[count++]);
            ++scaledIt;
        }

        // SamplEn2D
        scaledSampEn2D->SetD(s);
        scaledSampEn2D->SetImage(scaledImage);
        scaledSampEn2D->ComputeEntropy();
        m_Entropy.push_back(scaledSampEn2D->GetEntropy());
    }
}


template< typename TInputImage >
void
ModifiedMultiscaleEntropy2DImageCalculator< TInputImage >
::SetRegion(const RegionType & region)
{
    m_Region = region;
    m_RegionSetByUser = true;
}


template< typename TInputImage >
void 
ModifiedMultiscaleEntropy2DImageCalculator< TInputImage >
::ParametersCertification()
{
    //    TODO Testar cada metodo de erro...nao aparece na saida e nao termina o codigo.
    if (InputImageDimension != 2) {
        itkWarningMacro( << "Wrong input image dimension: "
                         << InputImageDimension << std::endl
                         << "This ITK class only accepts 2D images as input data.");
    }

    if ( m_M < 1 )
    {
        itkWarningMacro( << "Wrong value for M parameters: "
                         << m_M << std::endl
                         << "M must be an integer positive value.");
    }
    if (m_R < 0)
    {
        itkWarningMacro( << "Wrong value for R parameters: "
                         << m_R << std::endl
                         << "R must be a float point positive value.");
    }

    if (m_S < 1)
    {
        itkWarningMacro( << "Wrong value for S parameters: "
                         << m_S << std::endl
                         << "S must be an integer positive value.");
    }

    //TODO -- B must be double
    // if (m_BGV < 0)
    // {
    //     itkWarningMacro( << "Wrong value for B parameters: "
    //                      << m_BGV << std::endl
    //                      << "B must be a float point positive value.");
    // }
}

template< typename TInputImage >
void
ModifiedMultiscaleEntropy2DImageCalculator< TInputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);

    os << indent << "SampEn2D: "
       << m_Entropy.data()
       << std::endl;
    itkPrintSelfObjectMacro( Image );
    os << indent << "Region: " << std::endl;
    m_Region.Print( os, indent.GetNextIndent() );
    os << indent << "Region set by User: " << m_RegionSetByUser << std::endl;
}
} // end namespace itk

#endif
