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
#include "itkShapedNeighborhoodIterator.h"
#include "itkStatisticsImageFilter.h"
#include "itkNumericTraits.h"


#include "stdlib.h"

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

template<typename TInputImage>
void SampEn2DImageCalculator<TInputImage>
::ActivateOffsetsShapedNeighborhoods(ShapeIteratorType *it, MParameterValueType mValue){
    for (int x = 0; x < mValue; ++x) {
        for (int y = 0; y < mValue; ++y) {
            it->ActivateOffset({{x, y}});
        }
    }
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
    radius[0] = static_cast<int>(m_M);
    radius[1] = static_cast<int>(m_M);
    ShapeIteratorType firstWindowItM(radius, m_Image, m_Region);
    ShapeIteratorType firstWindowItMp1(radius, m_Image, m_Region);
    ShapeIteratorType secondWindowItM(radius, m_Image, m_Region);
    ShapeIteratorType secondWindowItMp1(radius, m_Image, m_Region);

    ActivateOffsetsShapedNeighborhoods(&firstWindowItM, m_M);
    ActivateOffsetsShapedNeighborhoods(&secondWindowItM,m_M);
    ActivateOffsetsShapedNeighborhoods(&firstWindowItMp1,m_M+1);
    ActivateOffsetsShapedNeighborhoods(&secondWindowItMp1,m_M+1);

    m_Entropy = NumericTraits< PixelType >::ZeroValue();

    typedef itk::StatisticsImageFilter<ImageType> StatisticsImageFilterType;
    typename StatisticsImageFilterType::Pointer statisticsImageFilter = StatisticsImageFilterType::New();
    statisticsImageFilter->SetInput(m_Image);
    statisticsImageFilter->Update();


    DoublePixelType SD = statisticsImageFilter->GetSigma();
    DoublePixelType tolerance = m_R*SD;

  //  cout<<"M: "<<m_M<<endl;
  //  cout<<"R: "<<m_R<<endl;
  //  cout<<"SD: "<<SD<<endl;
  //  cout<<"Tolerance: "<<tolerance<<endl;

    /** Matches for (m+1)-length patterns */
    MParameterValueType A = NumericTraits<MParameterValueType>::ZeroValue();
    /** Matches for (m)-length patterns */
    MParameterValueType B = NumericTraits<MParameterValueType>::ZeroValue();


    MParameterValueType Cim= NumericTraits<MParameterValueType>::ZeroValue();
    MParameterValueType Cim1= NumericTraits<MParameterValueType>::ZeroValue();
    RParameterValueType CimDensity= NumericTraits<RParameterValueType>::ZeroValue();
    RParameterValueType Cim1Density= NumericTraits<RParameterValueType>::ZeroValue();

    typename ImageType::RegionType region = m_Image->GetLargestPossibleRegion();
    typename ImageType::SizeType size = region.GetSize();
    RParameterValueType density = size[0]*size[1];

    firstWindowItM.Begin();
    firstWindowItMp1.Begin();
    while(!firstWindowItM.IsAtEnd()){
        secondWindowItM.Begin();
        secondWindowItMp1.Begin();
        while(!secondWindowItM.IsAtEnd()){
            //            firstWindowIt.Get
            if(IsSimilar(&firstWindowItM, &secondWindowItM,tolerance)){
                ++B;
                ++Cim;
            }
            if(IsSimilar(&firstWindowItMp1, &secondWindowItMp1,tolerance)){
                ++A;
                ++Cim1;
            }
            ++secondWindowItM;
            ++secondWindowItMp1;
        }
        CimDensity += Cim / (density - 1);
        Cim1Density += Cim1 / (density - 1);
        ++firstWindowItM;
        ++firstWindowItMp1;
    }

 //   cout<<"density: "<<density<<endl;

    CimDensity /= density;
    Cim1Density /= density;

  //  cout<<"CimDensity: "<<CimDensity<<endl;
  //  cout<<"Cim1Density: "<<Cim1Density<<endl;

    m_Entropy = static_cast<double>(-1)*std::log(static_cast<double>(Cim1Density)/static_cast<double>(CimDensity));
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
::IsSimilar(ShapeIteratorType *it1, ShapeIteratorType *it2, DoublePixelType tolerance)
{
    double diff = static_cast<double>(0);

    typename ShapeIteratorType::ConstIterator innerIt1 = it1->Begin();
    typename ShapeIteratorType::ConstIterator innerIt2 = it2->Begin();

    while (!innerIt1.IsAtEnd()) {
        diff = innerIt1.Get()-innerIt2.Get();
        if(diff >= tolerance){
            return false;
        }
        ++innerIt1;
        ++innerIt2;
    }

    return true;
}

template< typename TInputImage >
void SampEn2DImageCalculator< TInputImage >
::ParametersCertification()
{
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
