#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkSampEn2DImageCalculator.h"

#include "stdlib.h"

using namespace std;
int main(int argc, char *argv[])
{
 if ( argc < 4 )
        {
          std::cerr << "Missing parameters. " << std::endl;
          std::cerr << "Usage: " << std::endl;
          std::cerr << argv[0]
                    << " inputImageFile M R"
                    << std::endl;
          return -1;
        }

    const unsigned int Dimension = 2;

    typedef unsigned char                       InputPixelType;
    typedef unsigned char                       OutputPixelType;
    typedef itk::Image<InputPixelType, Dimension>    InputImageType;

    typedef itk::ImageFileReader<InputImageType> ReaderType;

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

    typedef itk::SampEn2DImageCalculator<InputImageType> SampEntropyType;
    SampEntropyType::Pointer sampEn2D = SampEntropyType::New();

    sampEn2D->SetImage(reader->GetOutput());
    sampEn2D->SetM(atoi(argv[2]));
    sampEn2D->SetR(atof(argv[3]));
    sampEn2D->ComputeEntropy();

    cout<<"SampEn2D: "<<sampEn2D->GetEntropy()<<endl;

  return EXIT_SUCCESS;
}
