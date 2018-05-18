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
    m_Entropy.clear();
    m_M = 1;
    m_R = 0.10;
    m_D = 1;
    m_BGV = 0.0;
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

    //Calculating stantard deviation only in the foreground area
    unsigned long long N = 0;
    double sigma = 0.0, mean = 0.0;
    ConstRegionIteratorType    copyIt(m_Image, m_Region);

    copyIt.GoToBegin(); // mean
    int count=0;
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
    sigma/=(N-1);
    sigma=sqrt(sigma);
    //End standard deviation

    //Assigning tolerance factor
    double tolerance = m_R*sigma;

    // Limit indexes
    size_t xLim = m_Nx - m_M*m_D; //col
    size_t yLim = m_Ny - m_M*m_D; //row

    //Counters of similar patterns for m and m+1 respectively
    unsigned long long B=0, A=0;
    
    // Total number of comparisons (paired patterns)
    N=0; 
    
    // Patterns reference indexes
    size_t p1, p2; 

    // Limits for patterns reference indexes
    size_t p1lim, p2lim;
    
    // Column and row iterators for template and moving patterns
    size_t i1, i2, j1, j2;

    //Pre-allocating index shifts for patterns matching
    int n_m = m_M*m_M;        // number of pixels of m length pattern
    int n_m1 = 2*m_M+1;     // number of extra pixels to complete m+1 length pattern
    int n_mm1 = n_m+n_m1; // number of pixels of m+1 length pattern
    
    size_t *s_m, *s_m1; //arrays with shifts to be applied for m and m+1 matching
    size_t *s_mm1; //array with shifts to be applied for background checking

    if ((s_m = (size_t *) calloc(n_m, sizeof(size_t))) == NULL) exit(1);
    if ((s_m1 = (size_t *) calloc(n_m1, sizeof(size_t))) == NULL) exit(1);  
    if ((s_mm1 = (size_t *) calloc(n_mm1, sizeof(size_t))) == NULL) exit(1);    
    
    int i,j; // aux indexes
    int kk;
    int k=0;
    for (j = 0; j < m_M ; j++){ // rows up to m
        kk = m_Nx * j * m_D;
        for (i = 0; i < m_M ; i++ ){ // cols up to m
            s_mm1[k] = s_m[k] = kk + i*m_D;
            k++;
        }
    }
    k=0;
    for (i = 0; i <= m_M; i++) {  // fill row m+1      
        s_mm1[n_m+k] = s_m1[k] = (m_Nx*m_M +i) * m_D;
        k++;
    }
    for (j = 0; j < m_M; j++) {  // fill col m+1
        s_mm1[n_m+k] = s_m1[k] = (m_Nx*j + m_M) * m_D;
        k++;
    }
    // End pre-allocating

    // Difference on pattern matching
    double diff;

    // Aux int variable
    int daux;


    // Allocating array to keep is-background info for all patterns
    bool NV[m_Nx*m_Ny];

    // Checking background
    for (j1 = 0; j1 < yLim; j1++) { // for each row
        p1lim = (j1*m_Nx)+m_Nx-(m_M*m_D); //set limit index for pattern index
        
        p1 = j1*m_Nx; //p1 is the pattern index
        
        while (p1 < p1lim) {
            // Check if pattern has any bg pixel
            NV[p1]=false;
            for (k=0; k < n_mm1; k++) {
                if(image_matrix[p1+s_mm1[k]] <= m_BGV) {
                    NV[p1]=true;
                    break; 
                }
            }
            p1++;
        }
    }
    // End background checking
        

    /* Starts running */
    for (j1 = 0; j1 < yLim; j1++) { // for each row
        p1lim = (j1*m_Nx)+m_Nx-(m_M*m_D); //set limit index for template pattern
        
        p1 = j1*m_Nx; //p1 is the template pattern
        i1 = 0;
        
        while (p1 < p1lim) {

            // if p1 is valid
            if(!NV[p1]) {

                //p2 is the moving pattern
                p2 = p1+m_D; //constraint to avoid spatial correlation effect
                i2 = i1+m_D;

                for (j2 = j1; j2 < yLim; j2++) { // for each row
                    p2lim = (j2*m_Nx)+m_Nx-(m_M*m_D); //set limit index for moving pattern
                                
                    while (p2 < p2lim ) {

                        // if p2 is valid
                        if(!NV[p2]) {

                            // Constraint to avoid spatial correlation effect
                            if ((j2-j1) < m_D) {
                                daux = (i1 > i2) ? i1 - i2 : i2 - i1;
                                if(daux < m_D) {
                                    i2++;  
                                    p2++;                  
                                    continue;
                                }
                            }

                            //increment number of paired patterns
                            N++;

                            //pattern matching for m with tolerance r
                            for (k=0; k < n_m; k++) {
                                diff = fabs(image_matrix[p1+s_m[k]] - image_matrix[p2+s_m[k]]);
                                if (diff >= tolerance) {
                                    diff = -1.0;
                                    break;
                                }
                            }

                            if (diff != -1.0) { // are patterns similar for m?
                                B++;
                                
                                //pattern matching for m+1 with tolerance r
                                for (k=0; k < n_m1; k++) {
                                    diff = fabs(image_matrix[p1+s_m1[k]] - image_matrix[p2+s_m1[k]]);
                                    if (diff >= tolerance) {
                                        diff = -1.0;
                                        break;
                                    }
                                }
                                
                                if (diff != -1.0) // are they similar for m+1?
                                    A++;
                            }// if m

                        }// if p2 valid

                        i2++;
                        p2++;
                    }//while p2

                    i2=0;
                    p2=(j2+1)*m_Nx;
                }//for j2

            }// if p1 valid

            i1++;
            p1++;
        } //while p1
    } //for j1

    /* Average probabilities Um and Um+1 may be obtained
       dividing B and A by N, respectively. */
    
    free(s_m);
    free(s_m1);
    free(s_mm1);

    m_Entropy.push_back(static_cast<double>(-1)*std::log(static_cast<double>(A)/static_cast<double>(B)));
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

    if (m_D < 1)
    {
        itkWarningMacro( << "Wrong value for D parameters: "
                         << m_D << std::endl
                         << "D must be an integer positive value.");
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
