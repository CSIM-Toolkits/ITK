#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkQualityImageFilter.h"
#include "itkCastImageFilter.h"

int main(int argc, char* argv[])
{
 if ( argc < 3 )
        {
          std::cerr << "Missing parameters. " << std::endl;
          std::cerr << "Usage: " << std::endl;
          std::cerr << argv[0]
                    << " inputReferenceImage inputCompareImage"
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
    typedef itk::ImageFileWriter<CastImageType> WriterType;

    ReaderType::Pointer readerReference = ReaderType::New();
    readerReference->SetFileName(argv[1]);
    ReaderType::Pointer readerCompare = ReaderType::New();
    readerCompare->SetFileName(argv[2]);


    try
        {
        readerReference->Update();
        readerCompare->Update();
        }
      catch ( itk::ExceptionObject &err)
        {
        std::cerr << "ExceptionObject caught !" << std::endl;
        std::cerr << err << std::endl;
        return -1;
        }

    typedef itk::QualityImageFilter<ImageType> FilterType;
    FilterType::Pointer filter = FilterType::New();
    filter->SetReferenceImage(readerReference->GetOutput());
    filter->SetCompareImage(readerCompare->GetOutput());

    std::cout<<filter->SNR()<<std::endl;
    std::cout<<filter->RMSE()<<std::endl;
    std::cout<<filter->SSIM()<<std::endl;


    return EXIT_SUCCESS;
}
