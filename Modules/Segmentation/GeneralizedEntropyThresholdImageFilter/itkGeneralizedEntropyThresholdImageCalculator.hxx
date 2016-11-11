
#ifndef __itkGeneralizedEntropyThresholdImageCalculator_txx
#define __itkGeneralizedEntropyThresholdImageCalculator_txx

#include "itkGeneralizedEntropyThresholdImageCalculator.h"
#include "itkImageRegionConstIteratorWithIndex.h"
#include "itkMinimumMaximumImageCalculator.h"

#include "vnl/vnl_math.h"

namespace itk
{ 
    
/**
 * Constructor
 */
template<class TInputImage>
GeneralizedEntropyThresholdImageCalculator<TInputImage>
::GeneralizedEntropyThresholdImageCalculator()
{
  m_Image = NULL;
  m_Threshold = NumericTraits<PixelType>::Zero;
  m_NumberOfHistogramBins = 128;
  m_RegionSetByUser = false;
  m_Q=1.0;
}


/*
 * Compute the MaxEntropy's threshold
 */
template<class TInputImage>
void
GeneralizedEntropyThresholdImageCalculator<TInputImage>
::Compute(void)
{
  if ( !m_Image ) { return; }
  if( !m_RegionSetByUser )
    {
    m_Region = m_Image->GetRequestedRegion();
    }

  double totalPixels = (double) m_Region.GetNumberOfPixels();
  if ( totalPixels == 0 ) { return; }


  // compute image max and min
  typedef MinimumMaximumImageCalculator<TInputImage> RangeCalculator;
  typename RangeCalculator::Pointer rangeCalculator = RangeCalculator::New();
  rangeCalculator->SetImage( m_Image );
  rangeCalculator->Compute();

  PixelType imageMin = rangeCalculator->GetMinimum();
  PixelType imageMax = rangeCalculator->GetMaximum();

  if ( imageMin >= imageMax )
    {
    m_Threshold = imageMin;
    return;
    }

  // create a histogram
  std::vector<double> relativeFrequency;
  relativeFrequency.resize( m_NumberOfHistogramBins );
  std::fill(relativeFrequency.begin(), relativeFrequency.end(), 0.0);

  double binMultiplier = (double) m_NumberOfHistogramBins /
    (double) ( imageMax - imageMin );

  typedef ImageRegionConstIteratorWithIndex<TInputImage> Iterator;
  Iterator iter( m_Image, m_Region );

  while ( !iter.IsAtEnd() )
    {
    unsigned int binNumber;
    PixelType value = iter.Get();

    if ( value == imageMin ) 
      {
      binNumber = 0;
      }
    else
      {
      binNumber = (unsigned int) vcl_ceil((value - imageMin) * binMultiplier ) - 1;
      if ( binNumber == m_NumberOfHistogramBins ) // in case of rounding errors
        {
        binNumber -= 1;
        }
      }

    relativeFrequency[binNumber] += 1.0;
    ++iter;

    }

  int threshold=-1;
  int ih, it;
  int first_bin;
  int last_bin;
  double tot_ent;  /* total entropy */
  double max_ent;  /* max entropy */
  double ent_back; /* entropy of the background pixels at a given threshold */
  double ent_obj;  /* entropy of the object pixels at a given threshold */
  std::vector<double> norm_histo(relativeFrequency.size()); /* normalized histogram */
  std::vector<double> P1(relativeFrequency.size()); /* cumulative normalized histogram */
  std::vector<double> P2(relativeFrequency.size());
  
  const double tolerance = 2.220446049250313E-16; // should get this
						  // from traits

  int total =0;
  for (ih = 0; (unsigned)ih < relativeFrequency.size(); ih++ )
    total+=relativeFrequency[ih];
  
  for (ih = 0; (unsigned)ih < relativeFrequency.size(); ih++ )
    norm_histo[ih] = (double)relativeFrequency[ih]/total;
  
  P1[0]=norm_histo[0];
  P2[0]=1.0-P1[0];
  for (ih = 1; (unsigned)ih < relativeFrequency.size(); ih++ )
    {
    P1[ih]= P1[ih-1] + norm_histo[ih];
    P2[ih]= 1.0 - P1[ih];
    }
  
  /* Determine the first non-zero bin */
  first_bin=0;
  for (ih = 0; (unsigned)ih < relativeFrequency.size(); ih++ ) 
    {
    if ( !(vcl_abs(P1[ih])<tolerance)) 
      {
      first_bin = ih;
      break;
      }
    }
  
  /* Determine the last non-zero bin */
  last_bin=relativeFrequency.size() - 1;
  for (ih = relativeFrequency.size() - 1; ih >= first_bin; ih-- ) 
    {
    if ( !(vcl_abs(P2[ih])<tolerance)) 
      {
      last_bin = ih;
      break;
      }
    }

  // Calculate the total entropy each gray-level
  // and find the threshold that maximizes it 
  max_ent = itk::NumericTraits<double>::min();
  
  for ( it = first_bin; it <= last_bin; it++ ) 
    {
    /* Entropy of the background pixels */
    ent_back = 0.0;
    for ( ih = 0; ih <= it; ih++ )  
      {
      if ( relativeFrequency[ih] !=0 ) 
	{
          if (m_Q==1.0) {
              ent_back -= ( norm_histo[ih] / P1[it] ) * vcl_log ( norm_histo[ih] / P1[it] );
          }else{
              ent_back += vcl_pow ( norm_histo[ih] / P1[it] , m_Q);
          }
	}
      }

  /* Entropy of the object pixels */
    ent_obj = 0.0;
    for ( ih = it + 1; (unsigned)ih < relativeFrequency.size(); ih++ )
      {
      if (relativeFrequency[ih]!=0)
	{
          if (m_Q==1.0) {
              ent_obj -= ( norm_histo[ih] / P2[it] ) * vcl_log ( norm_histo[ih] / P2[it] );
          }else{
              ent_obj += vcl_pow ( norm_histo[ih] / P2[it] , m_Q);
          }
	}
      }
    
  /* Total entropy */
    if (m_Q==1.0) {
        tot_ent = ent_back + ent_obj;
    }else{
        tot_ent = (1.0-ent_back)/( m_Q - 1.0 ) + (1.0-ent_obj)/( m_Q - 1.0 );
    }
  
  // IJ.log(""+max_ent+"  "+tot_ent);
  if ( max_ent < tot_ent ) 
    {
    max_ent = tot_ent;
    threshold = it;
    }
    }
  
  
  m_Threshold = static_cast<PixelType>( imageMin + 
                                        ( threshold ) / binMultiplier );


}

template<class TInputImage>
void
GeneralizedEntropyThresholdImageCalculator<TInputImage>
::SetRegion( const RegionType & region )
{
  m_Region = region;
  m_RegionSetByUser = true;
}

  
template<class TInputImage>
void
GeneralizedEntropyThresholdImageCalculator<TInputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf(os,indent);

  os << indent << "Q: " << m_Q << std::endl;
  os << indent << "Threshold: " << m_Threshold << std::endl;
  os << indent << "NumberOfHistogramBins: " << m_NumberOfHistogramBins << std::endl;
  os << indent << "Image: " << m_Image.GetPointer() << std::endl;
}

} // end namespace itk

#endif
