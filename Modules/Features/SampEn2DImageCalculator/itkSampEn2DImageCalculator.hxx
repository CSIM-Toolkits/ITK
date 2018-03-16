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
    m_Entropy = 0.0;
    m_M = 1;
    m_R = 0.10;
    m_D = 1;
    m_BackgroundValue = 0.0;
    m_RegionSetByUser = false;
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

    typename ImageType::SizeType size = m_Region.GetSize();

    m_Nx = size[0]; //ncols
    m_Ny = size[1]; //nrows

    double image_matrix[m_Nx*m_Ny];
    double sigma = 0.0, mean = 0.0;
    unsigned int N = 0;

    ConstRegionIteratorType    copyIt(m_Image, m_Region);
    copyIt.GoToBegin();
    int count=0;
    while (!copyIt.IsAtEnd()) {
        image_matrix[count++]=copyIt.Get();

        if (copyIt.Get()!=m_BackgroundValue) {
            N++;
            mean+=copyIt.Get();
        }
        ++copyIt;
    }
    mean/=N;

    //Calculating stantard deviation only in the foreground area
    copyIt.GoToBegin();
    while (!copyIt.IsAtEnd()) {
        if (copyIt.Get()!=m_BackgroundValue) {
            sigma+=pow(copyIt.Get()-mean,2.0);
        }
        ++copyIt;
    }
    sigma/=(N-1);
    sigma=sqrt(sigma);

    DoublePixelType tolerance = m_R*sigma;

    //Limit indexes for cols and rows
    int xLim = m_Nx - m_M*m_D; 
    int yLim = m_Ny - m_M*m_D; 
    
    // Total number of valid patterns (for both m and m+1)
    double den=0; 
   
    //Counters of similar patterns for m and m+1
    unsigned int B, A;  // for each template pattern   
    unsigned long long Cim=0, Cim1=0;  //sum over all template patters

    //Probabilities of patterns ocurrence for m and m+1
    double Cm, Cm1; 
    
    
    int i, j; //iterators 
    int ilim, jlim; // limits for pattern's indexes
        
    //Pre-allocating indexes shift for patterns matching and validation
    int nm = m_M*m_M;  //number of pixels in pattern m
    int nm1 = 2*m_M+1; //additional number of pixels for m+1 given m
    int nmm1 = nm+nm1; //number of pixels in pattern m+1
    
    int *kdm, *kdm1, *kdmm1; //arrays with shift indexes for m and m+1 matching and validation
    if ((kdm = (int *) calloc(nm, sizeof(int))) == NULL) exit(1);
    if ((kdm1 = (int *) calloc(nm1, sizeof(int))) == NULL) exit(1);
    if ((kdmm1 = (int *) calloc(nmm1, sizeof(int))) == NULL) exit(1);
    
    int kk;  
    int k=0; 
    for (i = 0; i < m_M ; i++){
        kk = m_Ny * i * m_D;
        for (j = 0; j < m_M ; j++ ){
            kdm[k] = kk + j*m_D;
            kdmm1[k] = kdm[k];
            k++;
        }
    }
    k=0;
    for (j = 0; j <= m_M; j++) {  // Compares column m+1      
        kdm1[k] = (m_Ny*m_M +j) * m_D;
        kdmm1[k+nm] = kdm1[k];
        k++;
    }
    for (i = 0; i < m_M; i++) {  // Compares row m+1
        kdm1[k] = (m_Ny*i + m_M) * m_D;
        kdmm1[k+nm] = kdm1[k];
        k++;
    }
    //End pre-allocating
    
    int ki, kj; //pattern's reference indexes
    double d; //differences during pattern matching
    bool valid; //valid pattern (has no background)
    double v;   //a moving pattern pixel value

    /* Starts running */
    for(i=0; i<yLim; i++) { //for each row
        ilim = (i*m_Nx)+xLim; //set limit index in current row for template pattern
        
        ki = i*m_Nx; //ki is top left pixel of template pattern
        
        while (ki < ilim) { //go across columns
            B = 0;
            A = 0;
          
            //check if it is a valid pattern
            valid = true;
            for (k=0; k < nmm1; k++) {
                if(image_matrix[ki+kdmm1[k]]==m_BackgroundValue) {
                    valid = false;
                    break; 
                }
            }
            
            if(valid) {//if it is a valid pattern (has no background)

                kj = ki+1; //kj is top left pixel of moving pattern

                for(j=i; j<yLim; j++) { //for each row
                    jlim = (j*m_Nx)+xLim; //set limit index in current row for moving pattern

                    while (kj < jlim ) { //go across columns
                        //pattern matching for m using tolerance
                        for (k=0; k < nm; k++) {
                            v = image_matrix[kj+kdm[k]];
                            if(v != m_BackgroundValue) {//if it is not a background pixel
                                d = fabs(image_matrix[ki+kdm[k]] - v);
                                if (d >= tolerance) {
                                    d = -1.0; //it is not a similar pattern
                                    break;
                                }
                            } else {
                                d = -1.0; //it is not a valid pattern (has background)
                                break;
                            }
                        }

                        // is pattern valid and similar for m?
                        if (d != -1.0) { 
                            valid=true;
                            //pattern matching for m+1 using tolerance
                            for (kk=0; kk < nm1; kk++) {
                                v = image_matrix[kj+kdm1[kk]];
                                if(v != m_BackgroundValue) {//if it is not a background pixel
                                    d = fabs(image_matrix[ki+kdm1[kk]] - v);
                                    if (d >= tolerance) {
                                        d = -1.0; //it is not a similar pattern
                                        break;
                                    }
                                } else {
                                    valid = false; //it is not a valid pattern (has background)
                                    break;
                                }
                            }

                            if (valid) { // is pattern valid so far?
                                //check remaining pixels
                                while(++kk < nm1) { 
                                    if(image_matrix[kj+kdm1[kk]] == m_BackgroundValue) //if it is a background pixel
                                      valid = false; //it is not a valid pattern
                                }
                                
                                // is pattern valid?
                                if (valid) {
                                    B++; //increment m counter
                                    if (d != -1.0) // is pattern similar for m+1?
                                        A++; //increment m+1 counter
                                }
                            }
                        }// if valid and similar for m

                        kj++;
                    }//while kj

                    kj= (j+1)*m_Nx;
                }//for j

                Cim += B;
                Cim1 += A;

                den++;
            }// if valid pattern
            
            ki++;
        } //while ki            
    } //for i
     
    /* outputs: probabilities */
    Cm = ((static_cast<double>(Cim)/(den-1.0))/den)*2.0;
    Cm1 = ((static_cast<double>(Cim1)/(den-1.0))/den)*2.0;
    
    free(kdm);
    free(kdm1);
    free(kdmm1);

    m_Entropy = static_cast<double>(-1)*std::log(static_cast<double>(Cm1)/static_cast<double>(Cm));
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
void 
SampEn2DImageCalculator< TInputImage >
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

    if (m_D < 1)
    {
        itkWarningMacro( << "Wrong value for D parameters: "
                         << m_D << std::endl
                         << "D must be an integer positive value.");
    }

//TODO -- B must be double
    // if (m_BackgroundValue < 0)
    // {
    //     itkWarningMacro( << "Wrong value for B parameters: "
    //                      << m_BackgroundValue << std::endl
    //                      << "B must be a float point positive value.");
    // }
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
