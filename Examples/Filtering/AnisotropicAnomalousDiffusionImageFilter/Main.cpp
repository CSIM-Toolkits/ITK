#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkAnisotropicAnomalousDiffusionImageFilter.h"
#include "itkCastImageFilter.h"

int main(int argc, char* argv[])
{
 if ( argc < 5 )
        {
          std::cerr << "Missing parameters. " << std::endl;
          std::cerr << "Usage: " << std::endl;
          std::cerr << argv[0]
                    << " inputImageFile outputImageFile Condutance Qvalue NumberOfIteration TimeStep"
                    << std::endl;
          return -1;
        }

    const unsigned int Dimension = 3;

    typedef unsigned char                       PixelType;
    typedef unsigned char                       PixelOutType;
    typedef itk::Image<PixelType, Dimension>    ImageType;
    typedef itk::Image<PixelOutType, Dimension> ImageOutType;
    typedef float                               CastPixelType;
    typedef itk::Image<CastPixelType, Dimension> CastImageType;

    typedef itk::ImageFileReader<ImageType> ReaderType;
    typedef itk::ImageFileWriter<ImageOutType> WriterType;

    ReaderType::Pointer reader = ReaderType::New();
    reader->SetFileName(argv[1]);
    try
        {
        reader->Update();
        }
      catch ( itk::ExceptionObject &err)
        {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return -1;
        }

    typedef itk::CastImageFilter<ImageType, CastImageType> CastType;
    typename CastType::Pointer cast = CastType::New();
    cast->SetInput(reader->GetOutput());
    cast->Update();

    typedef itk::AnisotropicAnomalousDiffusionImageFilter<CastImageType, CastImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();

    filter->SetInput(cast->GetOutput());
    filter->SetCondutance(std::atof(argv[3]));
    filter->SetQ(std::atof(argv[4]));
    filter->SetIterations(std::atoi(argv[5]));
    filter->SetTimeStep(std::atof(argv[6]));
    filter->Update();

    typedef itk::CastImageFilter<CastImageType, ImageOutType> CastBackType;
    typename CastBackType::Pointer castBack = CastBackType::New();
    castBack->SetInput(filter->GetOutput());
    castBack->Update();


      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName(argv[2]);
      writer->SetInput( castBack->GetOutput() );


      try
        {
        writer->Update();
        }
      catch ( itk::ExceptionObject &err)
        {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return -1;
        }

  return EXIT_SUCCESS;
}
