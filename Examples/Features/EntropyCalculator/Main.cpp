#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"

#include "itkModifiedMultiscaleEntropy2DImageCalculator.h"
#include "itkSampEn2DImageCalculator.h"

#include "stdlib.h"

using namespace std;
int main(int argc, char *argv[])
{
 if ( argc < 6 )
        {
          std::cerr << "Missing parameters. " << std::endl;
          std::cerr << "Usage: " << std::endl;
          std::cerr << argv[0]
                    << " inputImageFile M R D BGV"
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

    typedef itk::SampEn2DImageCalculator<InputImageType> SampEn2DType;
    SampEn2DType::Pointer sampEn2D = SampEn2DType::New();
    sampEn2D->UseRParameterAsPercentageOn();
    sampEn2D->SetImage(reader->GetOutput());
    sampEn2D->SetM(atoi(argv[2]));
    sampEn2D->SetR(atof(argv[3]));
    sampEn2D->SetD(atoi(argv[4]));
    sampEn2D->SetBGV(atof(argv[5]));
    sampEn2D->ComputeEntropy();

    cout<<"SampEn2D: "<<sampEn2D->GetEntropy()<<endl;

    typedef itk::ModifiedMultiscaleEntropy2DImageCalculator<InputImageType> MMSE2DType;
    MMSE2DType::Pointer mmse2D = MMSE2DType::New();

    mmse2D->SetImage(reader->GetOutput());
    mmse2D->SetM(atoi(argv[2]));
    mmse2D->SetR(atof(argv[3]));
    mmse2D->SetD(atoi(argv[4]));
    mmse2D->SetBGV(atof(argv[5]));
    mmse2D->ComputeEntropy();

    cout<<"MMSE2D: "<<mmse2D->GetEntropy().data()<<endl;

  return EXIT_SUCCESS;
}
