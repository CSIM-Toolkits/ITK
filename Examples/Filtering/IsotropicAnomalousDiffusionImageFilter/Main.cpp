#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "IsotropicAnomalousDiffusionImageFilter.h"

int main(int, char*[])
{
  // Setup types
  typedef itk::Image<unsigned char, 2>   ImageType;
  typedef itk::ImageFileReader<ImageType> ReaderType;
  typedef itk::ImageFileWriter<ImageType> WriterType;
    typedef itk::IsotropicAnomalousDiffusionImageFilter<ImageType, ImageType>  FilterType;

  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName("bizu.png");//Set other image filename. It's better if the image file is already float pixel data. 
  reader->Update();

  // Create and set the filter parameters
  FilterType::Pointer filter = FilterType::New();
  filter->SetInput(reader->GetOutput());
  filter->SetIterations(1);
  filter->SetQ(1.0);
  filter->Update();

  WriterType::Pointer writer = WriterType::New();
  writer->SetFileName("bizu-out.tif");//Set an output filename. Tif format file can save float pixel data.
  writer->SetInput(filter->GetOutput());
  writer->Update();

  return EXIT_SUCCESS;
}
