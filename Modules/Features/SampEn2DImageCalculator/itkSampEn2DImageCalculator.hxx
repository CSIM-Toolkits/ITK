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
#ifndef __itkSampEn2DImageCalculator_hxx
#define __itkSampEn2DImageCalculator_hxx

#include "itkSampEn2DImageCalculator.h"

using namespace std;

namespace itk
{
/**
 * Constructor
 */
template< typename TInputImage >
SampEn2DImageCalculator< TInputImage >
::SampEn2DImageCalculator()
{
    m_Image = TInputImage::New();
    m_Entropy = NumericTraits< PixelType >::ZeroValue();
    m_M = 1;
    m_R = 0.10;
    m_RegionSetByUser = false;
    m_IsSimilar = false;
}

/**
 * Compute Entropy of m_Image
 */
template< typename TInputImage >
void
SampEn2DImageCalculator< TInputImage >
::ComputeEntropy(void)
{
    if ( !m_RegionSetByUser )
    {
        m_Region = m_Image->GetRequestedRegion();
    }

    ParametersCertification();

    typename ImageType::SizeType radius;
    radius.Fill(static_cast<int>(m_M));

    m_Entropy = NumericTraits< PixelType >::ZeroValue();

    typedef itk::StatisticsImageFilter<ImageType> StatisticsImageFilterType;
    typename StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
    statisticsImageFilter->SetInput(m_Image);
    statisticsImageFilter->Update();

    DoublePixelType tolerance = m_R*statisticsImageFilter->GetSigma();

    typename ImageType::SizeType size = m_Region.GetSize();

    //    clock_t begin, end;
    //    begin=clock();

    m_Nx = size[0];
    m_Ny = size[1];

    double image_matrix[m_Nx*m_Ny];

    ConstRegionIteratorType    copyIt(m_Image, m_Region);
    copyIt.GoToBegin();
    int count=0;
    while (!copyIt.IsAtEnd()) {
        image_matrix[count++]=copyIt.Get();
        ++copyIt;
    }


    int Cim, Cim1;
    double Cm = 0.0;
    double Cm1 = 0.0;
    // Total number of patterns (for both m and m+1)
    double den = (m_Nx - m_M) * (m_Ny - m_M);

    for (int yi = 0; yi < m_Ny - m_M; yi++) {
        for (int xi = 0; xi < m_Nx - m_M; xi++) {

            // Counters of similar patterns for m and m+1
            Cim = Cim1 = 0;

            int yj = yi;
            int xj = xi + 1;
            while (xj < m_Nx - m_M) {
                if (similar(image_matrix, xi, yi, xj, yj, m_M, tolerance)) {  // Similar for M?
                    Cim++;

                    // Are they still similar for the next point?
                    if (similarNext(image_matrix, xi, yi, xj, yj, m_M, tolerance)) { // Similar for M?
                        Cim1++;
                    }
                }
                xj++;
            }

            for (yj = yi + 1; yj < m_Ny - m_M; yj++) {
                for (xj = 0; xj < m_Nx - m_M; xj++) {
                    if (similar(image_matrix, xi, yi, xj, yj, m_M, tolerance)) {  // Similar for M?
                        Cim++;

                        // Are they still similar for the next point?
                        if (similarNext(image_matrix, xi, yi, xj, yj, m_M, tolerance)) { // Similar for M?
                            Cim1++;
                        }
                    }
                }
            }

            Cm += Cim / (den - 1);
            Cm1 += Cim1 / (den - 1);
        }
    }
    Cm /= den;
    Cm1 /= den;

    m_Entropy = static_cast<double>(-1)*std::log(static_cast<double>(Cm1)/static_cast<double>(Cm));

    //    end=clock();
    //    double time = (double)(end-begin)/CLOCKS_PER_SEC;
    //    cout<<"time "<<time<<endl;

}


template< typename TInputImage >
void
SampEn2DImageCalculator< TInputImage >
::SetRegion(const RegionType & region)
{
    m_Region = region;
    m_RegionSetByUser = true;
}

template< typename TInputImage >
bool
SampEn2DImageCalculator< TInputImage >
::similar(double* image, int x1, int y1, int x2, int y2, int m, double r)
{
    for (int y = 0; y < m; y++) {
        for (int x = 0; x < m; x++) {
            double diff = std::abs(image[x1 + x + (y1 + y)*m_Nx] - image[x2 + x + (y2 + y)*m_Nx]);
            if (diff >= r) {
                return false;
            }
        }
    }
    return true;
}

template< typename TInputImage >
bool
SampEn2DImageCalculator< TInputImage >
::similarNext(double* image, int x1, int y1, int x2, int y2, int m, double r)
{
    double diff;
    for (int y = 0; y <= m; y++) {  // Compares collumn M
        diff = std::abs(image[x1 + m + (y1 + y)*m_Nx] - image[x2 + m + (y2 + y)*m_Nx]);
        if (diff >= r) {
            return false;
        }
    }

    for (int x = 0; x <= m; x++) {  // Compares row M
        diff = std::abs(image[x1 + x + (y1 + m)*m_Nx] - image[x2 + x + (y2 + m)*m_Nx]);
        if (diff >= r) {
            return false;
        }
    }

    return true;
}

template< typename TInputImage >
void SampEn2DImageCalculator< TInputImage >
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
    if (m_R < 0){
        itkWarningMacro( << "Wrong value for R parameters: "
                         << m_R << std::endl
                         << "R must be a float point positive value.");
    }
}

template< typename TInputImage >
void
SampEn2DImageCalculator< TInputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
    Superclass::PrintSelf(os, indent);

    os << indent << "SampEn2D: "
       << static_cast< typename NumericTraits< PixelType >::PrintType >( m_Entropy )
       << std::endl;
    itkPrintSelfObjectMacro( Image );
    os << indent << "Region: " << std::endl;
    m_Region.Print( os, indent.GetNextIndent() );
    os << indent << "Region set by User: " << m_RegionSetByUser << std::endl;
}
} // end namespace itk

#endif
