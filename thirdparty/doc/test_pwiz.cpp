//#include "pwiz_tools/common/FullReaderList.hpp"
#include "pwiz/data/msdata/DefaultReaderList.hpp"
#include "pwiz/data/msdata/MSDataFile.hpp"
#include "pwiz/utility/misc/Std.hpp"
#include "pwiz/utility/misc/Filesystem.hpp"

using namespace pwiz::msdata;
using namespace pwiz::cv;

int main(int argc, char* argv[])
{
	if (argc == 1)
	{
		cerr << "Usage: pwiz_example <MS data file> <output .mzXML filename>" << endl;
		return 1;
	}

	try
	{
		//use both opern format readers and vendor readers;
		//on some platforms the vendor readers can identify vendor data but not read it
		
    DefaultReaderList readers;

		//populate an MSData object from an MS data filepath
		
    MSDataFile msd(argv[1], &readers);

		SpectrumList& spectrumList = *msd.run.spectrumListPtr;
		SpectrumPtr spectrum;
		//Read m/z and intensity values from the spectra
		const bool getBinaryData = true;	
		size_t numSpectra = spectrumList.size();
		if (numSpectra < 2){
			throw runtime_error("Input file must contain at least two spectra!");
		}
		for (int i = 0; i<2; ++i){
			spectrum = spectrumList.spectrum(i, getBinaryData);
			vector<MZIntensityPair> pairs;
			spectrum->getMZIntensityPairs(pairs);
			//now pairs is a list of all of the <m/z,intensity> peaks in the spectrum
			int counter = 0;
			//iterate through the first 5 peaks (at most) in the spectrum, and output them to stdout
			for (vector<MZIntensityPair>::const_iterator it = pairs.begin(), end = pairs.end(); it!=end; ++it)
			{
				if (counter > 4){
					break;
				}
				cout<<it->mz<<'\t'<< it->intensity<<endl;
				++counter;
			}
			cout<<"-----------------------------------"<<endl;
		}

		//Now, a demonstration on how to output a file, we'll output an mzML
		//MSDataFile::WriteConfig holds information on how the output file will be written
		MSDataFile::WriteConfig writeConfig;
		//we'd like to write an mzML file:
		writeConfig.format = MSDataFile::Format_mzML;
		//we can set the precision of the output data if we want, I'll write the defaults for demonstration purposes
		writeConfig.binaryDataEncoderConfig.precision = BinaryDataEncoder::Precision_64;
		writeConfig.binaryDataEncoderConfig.precisionOverrides[MS_m_z_array] = BinaryDataEncoder::Precision_64;
		writeConfig.binaryDataEncoderConfig.precisionOverrides[MS_intensity_array] = BinaryDataEncoder::Precision_32;

		//write the file
		MSDataFile::write(msd, argv[2], writeConfig);
	}
	catch (runtime_error &e)
	{
		cerr <<"Caught exception: "<<e.what() <<endl;
		return 1;
	}
	return 0;
}
