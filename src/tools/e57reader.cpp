#include <array>
#include <fstream>
#include <experimental/filesystem>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <cmath>

#include "E57Foundation.h"
#include "E57Simple.h"

double dot4(const std::array<double, 4> &a, const std::array<double, 4> &b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3];
}

void printNames(int nbParam, char *param[]);
void getGroupsOBJ(int nbParam, char *param[]);
void e57_to_bin(int nbParam, char *param[]);
void oneStation_E57_to_bin(e57::Reader &eReader, int scanIndex, char *param[]);
template <class Matrix, class Points, class Colors>
void write_station_to_bin(const std::string &ptx_fn,
                          int nColumn, int nRow,
                          const Matrix &matrix,
                          const std::vector<Points> &points,
                          const std::vector<Colors> &colors,
                          const std::vector<double> &intensity);
template <class Matrix, class Points, class Colors>
void write_station_to_ptx(const std::string &ptx_fn,
                          int nColumn, int nRow,
                          const Matrix &matrix,
                          const std::vector<Points> &points,
                          const std::vector<Colors> &colors,
                          const std::vector<double> &intensity);

int main(int nbParam, char *param[])
{
    try
    {
        //printNames(nbParam, param);
        e57_to_bin(nbParam, param);
    }
    catch (e57::E57Exception& ex) {
        ex.report(__FILE__, __LINE__, __FUNCTION__);
    }
    catch (std::exception& ex) {
        cerr << "Got an std::exception, what=" << ex.what() << endl;
    }
    catch (...) {
        cerr << "Got an unknown exception" << endl;
    }

    return 0;
}

void e57_to_bin(int nbParam, char *param[])
{
        // Create a ReaderImpl

    e57::ustring bsFile(param[1]);			//converts Unicode to UTF-8
    e57::Reader	eReader(bsFile);

    ///////////////////////////////////////////////////////////////
    // ACCESSING ROOT DATA

    // Read the root

    e57::E57Root	rootHeader;
    eReader.GetE57Root(rootHeader);

    //Access all the root information like

    e57::ustring &fileGuid = rootHeader.guid;
    double fileGPSTime = rootHeader.creationDateTime.dateTimeValue;

    unsigned short	utc_year;	//Universal Time Coordinated [year]
    unsigned char	utc_month;	//1-12 months
    unsigned char	utc_day;	//1-31 days
    unsigned char	utc_hour;	//hours
    unsigned char	utc_minute;	//minutes
    float		utc_seconds;//seconds
    unsigned short	gps_week;	//GPS week (0-1024+)
    double		gps_tow;	//GPS time of week(0-604800.0) seconds

   // e57::GetUTCFromGPSDateTime(
   //     fileGPSTime,		//!< GPS Date Time
   //     &utc_year,		//!< The year 1900-9999
   //     &utc_month,		//!< The month 0-11
   //     &utc_day,		//!< The day 1-31
   //     &utc_hour,		//!< The hour 0-23
   //     &utc_minute,		//!< The minute 0-59
   //     &utc_seconds	//!< The seconds 0.0 - 59.999
   // )

    gps_week = ((int)floor(fileGPSTime)) / 604800;
    gps_tow = fileGPSTime - gps_week*604800.;

    bool result = TIMECONV_GetUTCTimeFromGPSTime(gps_week, gps_tow,
        &utc_year, &utc_month, &utc_day, &utc_hour, &utc_minute, &utc_seconds);

    ///////////////////////////////////////////////////////////////
    // ACCESSING SCAN DATA3D

    //Get the number of scan images available

    int data3DCount = eReader.GetData3DCount();
    for (int scanIndex = 0; scanIndex < data3DCount; ++scanIndex)
    {
        oneStation_E57_to_bin(eReader, scanIndex, param);
    }

    eReader.Close();
}

void oneStation_E57_to_bin(e57::Reader &eReader, int scanIndex, char *param[])
{

    //Read the scan 0 header.
    e57::Data3D		scanHeader;
    eReader.ReadData3D(scanIndex, scanHeader);
    //Access all the header information

    e57::ustring bsFile(param[1]);
    e57::ustring bstrName = scanHeader.name.c_str();
    e57::ustring bstrGuid = scanHeader.guid.c_str();
    e57::ustring bstrDesc = scanHeader.description.c_str();

    unsigned short	utc2_year;	//Universal Time Coordinated [year]
    unsigned char	utc2_month;	//1-12 months
    unsigned char	utc2_day;	//1-31 days
    unsigned char	utc2_hour;	//hours
    unsigned char	utc2_minute;	//minutes
    float		utc2_seconds;//seconds
    unsigned short	gps2_week;	//GPS week (0-1024+)
    double		gps2_tow;	//GPS time of week(0-604800.0) seconds

    double startGPSTime = scanHeader.acquisitionStart.dateTimeValue;
    gps2_week = ((int)floor(startGPSTime)) / 604800;
    gps2_tow = startGPSTime - gps2_week*604800.;

    bool result = TIMECONV_GetUTCTimeFromGPSTime(gps2_week, gps2_tow,
        &utc2_year, &utc2_month, &utc2_day, &utc2_hour, &utc2_minute, &utc2_seconds);
    double endGPSTime = scanHeader.acquisitionEnd.dateTimeValue;

    //Get pose information

    struct {
        float m00, m01, m02, m03;
        float m10, m11, m12, m13;
        float m20, m21, m22, m23;
        float m30, m31, m32, m33;
    } matrix;

    double X = scanHeader.pose.rotation.x,
        Y = scanHeader.pose.rotation.y,
        Z = scanHeader.pose.rotation.z,
        W = scanHeader.pose.rotation.w;
    double xx = X * X;
    double xy = X * Y;
    double xz = X * Z;
    double xw = X * W;

    double yy = Y * Y;
    double yz = Y * Z;
    double yw = Y * W;

    double zz = Z * Z;
    double zw = Z * W;

    matrix.m00 = 1 - 2 * (yy + zz);
    matrix.m01 = 2 * (xy - zw);
    matrix.m02 = 2 * (xz + yw);

    matrix.m10 = 2 * (xy + zw);
    matrix.m11 = 1 - 2 * (xx + zz);
    matrix.m12 = 2 * (yz - xw);

    matrix.m20 = 2 * (xz - yw);
    matrix.m21 = 2 * (yz + xw);
    matrix.m22 = 1 - 2 * (xx + yy);


    matrix.m03 = scanHeader.pose.translation.x;
    matrix.m13 = scanHeader.pose.translation.y;
    matrix.m23 = scanHeader.pose.translation.z;
    matrix.m30 = matrix.m31 = matrix.m32 = 0;
    matrix.m33 = 1;

    auto dirName = param[2] + scanHeader.name;///*bsFile.substr(0, 9) + "_" +*/ bstrDesc.substr(11, 11);//
    std::string defaultName = dirName + "/" + scanHeader.name;///*bsFile.substr(0, 9) + "_" +*/ bstrDesc.substr(11, 11);//
    std::string station_name = defaultName;
    {//check if a folder with the same name already exists. If so, append _bis-x, where x is a number from 1 up to find a different name
        auto dirNamebkp = dirName;
        auto station_namebkp = station_name;
        int count = 1;
        while (!std::experimental::filesystem::create_directory(dirName))
        {
            dirName = dirNamebkp + "_bis-" + std::to_string(count);
            station_name = station_namebkp + "_bis-" + std::to_string(count);
            ++count;
        }
    }
    
    std::cout << station_name.substr(station_name.find_last_of('/')+1, station_name.size());
    //   return 0;
    //   ISI::Point translation;
    //   translation.x(scanHeader.pose.translation.x);
    //   translation.y(scanHeader.pose.translation.y);
    //   translation.z(scanHeader.pose.translation.z);
    //
    //   e57::ISI::Quat rotation;
    //   rotation.w(scanHeader.pose.rotation.w);
    //   rotation.x(scanHeader.pose.rotation.x);
    //   rotation.y(scanHeader.pose.rotation.y);
    //   rotation.z(scanHeader.pose.rotation.z);

    //Get scanner information

    e57::ustring bstrSerial = scanHeader.sensorSerialNumber.c_str();
    e57::ustring bstrVendor = scanHeader.sensorVendor.c_str();
    e57::ustring bstrModel = scanHeader.sensorModel.c_str();
    e57::ustring bstrSoftware = scanHeader.sensorSoftwareVersion.c_str();
    e57::ustring bstrFirmware = scanHeader.sensorFirmwareVersion.c_str();
    e57::ustring bstrHardware = scanHeader.sensorHardwareVersion.c_str();

    //Get environment information

    double temperature = scanHeader.temperature;
    double humidity = scanHeader.relativeHumidity;
    double airPressure = scanHeader.atmosphericPressure;

    /////////////////////////////////////////////////////////////////
    // ACCESSING SCAN DATA

    //Get the Size of the Scan
    int64_t nColumn = 0;
    int64_t nRow = 0;

    int64_t nPointsSize = 0;	//Number of points

    int64_t nGroupsSize = 0;	//Number of groups
    int64_t nCountSize = 0;		//Number of points per group
    bool	bColumnIndex = false; //indicates that idElementName is "columnIndex"

    eReader.GetData3DSizes(scanIndex, nRow, nColumn, nPointsSize, nGroupsSize, nCountSize, bColumnIndex);
    std::cout << "\t\t" << nColumn << "x" << nRow << std::endl;
    int64_t nSize = nRow;
    if (nSize == 0) nSize = 1024;	// choose a chunk size

                                    //Setup buffers

    std::vector<int8_t> isInvalidData;
    if (scanHeader.pointFields.cartesianInvalidStateField)
        isInvalidData.resize(nSize);

    //Setup Points Buffers

    bool sphericalMode = false;
    //no cartesian fields?
    if (!scanHeader.pointFields.cartesianXField &&
        !scanHeader.pointFields.cartesianYField &&
        !scanHeader.pointFields.cartesianZField)
    {
        //let's look for spherical ones
        if (!scanHeader.pointFields.sphericalRangeField &&
            !scanHeader.pointFields.sphericalAzimuthField &&
            !scanHeader.pointFields.sphericalElevationField)
        {
            throw "No readable point in scan " + station_name;
        }
        sphericalMode = true;
    }

    std::vector<double> xData;
    if (scanHeader.pointFields.cartesianXField)
        xData.resize(nSize);

    std::vector<double> yData;
    if (scanHeader.pointFields.cartesianYField)
        yData.resize(nSize);

    std::vector<double> zData;
    if (scanHeader.pointFields.cartesianZField)
        zData.resize(nSize);

    std::vector<double> range;
    if (scanHeader.pointFields.sphericalRangeField)
        range.resize(nSize);

    std::vector<double> azimuth;
    if (scanHeader.pointFields.sphericalAzimuthField)
        azimuth.resize(nSize);

    std::vector<double> elevation;
    if (scanHeader.pointFields.sphericalElevationField)
        elevation.resize(nSize);

    //Setup intensity buffers if present

    std::vector<double> intData;
    bool		bIntensity = false;
    double		intRange = 0;
    double		intOffset = 0;

    if (scanHeader.pointFields.intensityField)
    {
        bIntensity = true;
        intData.resize(nSize);
        intRange = scanHeader.intensityLimits.intensityMaximum - scanHeader.intensityLimits.intensityMinimum;
        intOffset = scanHeader.intensityLimits.intensityMinimum;
    }

    //Setup color buffers if present
    
    std::vector<uint16_t> redData;
    std::vector<uint16_t> greenData;
    std::vector<uint16_t> blueData;
    bool		bColor = false;
    int32_t		colorRedRange = 1;
    int32_t		colorRedOffset = 0;
    int32_t		colorGreenRange = 1;
    int32_t		colorGreenOffset = 0;
    int32_t		colorBlueRange = 1;
    int32_t		colorBlueOffset = 0;

    if (scanHeader.pointFields.colorRedField)
    {
        bColor = true;
        redData.resize(nSize);
        greenData.resize(nSize);
        blueData.resize(nSize);
        colorRedRange = scanHeader.colorLimits.colorRedMaximum - scanHeader.colorLimits.colorRedMinimum;
        colorRedOffset = scanHeader.colorLimits.colorRedMinimum;
        colorGreenRange = scanHeader.colorLimits.colorGreenMaximum - scanHeader.colorLimits.colorGreenMinimum;
        colorGreenOffset = scanHeader.colorLimits.colorGreenMinimum;
        colorBlueRange = scanHeader.colorLimits.colorBlueMaximum - scanHeader.colorLimits.colorBlueMinimum;
        colorBlueOffset = scanHeader.colorLimits.colorBlueMinimum;
    }

    //Setup the GroupByLine buffers information

    std::vector<int64_t> idElementValue;
    std::vector<int64_t> startPointIndex;
    std::vector<int64_t> pointCount;
    if (nGroupsSize > 0)
    {
        idElementValue.resize(nGroupsSize);
        startPointIndex.resize(nGroupsSize);
        pointCount.resize(nGroupsSize);

        if (!eReader.ReadData3DGroupsData(scanIndex, nGroupsSize, idElementValue.data(),
            startPointIndex.data(), pointCount.data()))
            nGroupsSize = 0;
    }

    //Setup row/column index information

    std::vector<int32_t> rowIndex;
    std::vector<int32_t> columnIndex;
    if (scanHeader.pointFields.rowIndexField)
        rowIndex.resize(nSize);
    if (scanHeader.pointFields.columnIndexField)
        columnIndex.resize(nRow);

    //Get dataReader object
#define GET_PTR(vector) vector.empty() ? nullptr : vector.data()
    e57::CompressedVectorReader dataReader = eReader.SetUpData3DPointsData(
        scanIndex,			//!< data block index given by the NewData3D
        nRow,				//!< size of each of the buffers given
        GET_PTR(xData),				//!< pointer to a buffer with the x data
        GET_PTR(yData),				//!< pointer to a buffer with the y data
        GET_PTR(zData),				//!< pointer to a buffer with the z data
        GET_PTR(isInvalidData),		//!< pointer to a buffer with the valid indication
        GET_PTR(intData),			//!< pointer to a buffer with the lidar return intesity
        NULL,
        GET_PTR(redData),			//!< pointer to a buffer with the color red data
        GET_PTR(greenData),			//!< pointer to a buffer with the color green data
        GET_PTR(blueData),			//!< pointer to a buffer with the color blue data
        NULL,
        GET_PTR(range),//range
        GET_PTR(azimuth),//Azimuth
        GET_PTR(elevation),//Elevation
        NULL,
        GET_PTR(rowIndex),			//!< pointer to a buffer with the rowIndex
        GET_PTR(columnIndex)			//!< pointer to a buffer with the columnIndex
    );

    //Read the point data

    int64_t		count = 0;
    unsigned	size = 0;
    int			col = 0;
    int			row = 0;

    uint64_t nbP = nRow * nColumn;
    uint64_t h = nRow;

    struct vec3
    {
        float x, y, z;
    };
    std::vector<vec3> points(nbP);
    struct colors_type
    {
        uint8_t  r, g, b;
    };
    std::vector<colors_type> colors(nbP);
    std::vector<double> intensity;
    if (bIntensity)
    {
        intensity.resize(nbP);
    }

    //std::cout << std::endl;
    //std::cout << "Reading " << param[1] << "...";
    while (size = dataReader.read())
    {
        std::stringstream ss_buffer;
        //Use the data
        // std::cout << "size: " << size << " : " << cnt << std::endl;
        for (long i = 0; i < size; i++)
        {
            if (!columnIndex.empty())
                col = columnIndex[i];
            else
                col = 0;	//point cloud case

            if (!rowIndex.empty())
                row = rowIndex[i];
            else
                row = count;	//point cloud case
            int idx = col*h + row;
            if (sphericalMode)
            {
                points[idx] = { float(range[i] * cos(elevation[i]) * cos(azimuth[i])),
                                float(range[i] * cos(elevation[i]) * sin(azimuth[i])),
                                float(range[i] * sin(elevation[i])) };
            }
            else
            {
                if (!isInvalidData.empty() && isInvalidData[i] == 0)
                    points[idx] = { float(xData[i]), float(yData[i]), float(zData[i]) };
                //pScan->SetPoint(row, col, xData[i], yData[i], zData[i]);
                else
                    col = 0;
            }
            if (bIntensity) {		//Normalize intensity to 0 - 1.
                intensity[idx] = (intData[i] - intOffset) / intRange;
                //    pScan->SetIntensity(row, col, intensity);
                if (!bColor) {
                    int y = 255. * std::max(std::min(intensity[idx], 1.), 0.);
                    colors[idx] = { uint8_t(y), uint8_t(y), uint8_t(y) };
                }
            }

            if (bColor) {			//Normalize color to 0 - 255
                int red = ((redData[i] - colorRedOffset) * 255) / colorRedRange;
                int green = ((greenData[i] - colorGreenOffset) * 255) / colorBlueRange;
                int blue = ((blueData[i] - colorBlueOffset) * 255) / colorBlueRange;
                //    pScan->SetColor(row, col, red, green, blue);
                colors[idx] = { uint8_t(red), uint8_t(green), uint8_t(blue) };
            }
            count++;
        }

    }
    //Close and clean up
    dataReader.close();

    ///////////////////////////////////////////////////////////////////////
    // ACCESSING PICTURE IMAGE2D

    //Get the number of picture images available

    int image2DCount = eReader.GetImage2DCount();
    for (int imageIndex = 0; imageIndex < image2DCount; ++imageIndex)
    {
        //selecting the first picture image

        //Read the picture 0 header.

        e57::Image2D imageHeader;
        eReader.ReadImage2D(imageIndex, imageHeader);

        //Access all the header information

        e57::ustring bstrName = imageHeader.name.c_str();
        e57::ustring bstrGuid = imageHeader.guid.c_str();
        e57::ustring bstrDesc = imageHeader.description.c_str();
        double imageGPSTime = imageHeader.acquisitionDateTime.dateTimeValue;

        //Get pose information

        //      ISI::Point translation;
        //      translation.x(imageHeader.pose.translation.x);
        //      translation.y(imageHeader.pose.translation.y);
        //      translation.z(imageHeader.pose.translation.z);
        //
        //      ISI::Quat rotation;
        //      rotation.w(imageHeader.pose.rotation.w);
        //      rotation.x(imageHeader.pose.rotation.x);
        //      rotation.y(imageHeader.pose.rotation.y);
        //      rotation.z(imageHeader.pose.rotation.z);

        //Get camera information

        e57::ustring bstrSerial = imageHeader.sensorSerialNumber.c_str();
        e57::ustring bstrVendor = imageHeader.sensorVendor.c_str();
        e57::ustring bstrModel = imageHeader.sensorModel.c_str();


        ///////////////////////////////////////////////////////////////////////
        // ACCESSING PICTURE IMAGE

        //Get the Size of the Picture

        e57::Image2DProjection	imageProjection;	//like E57_SPHERICAL
        e57::Image2DType		imageType;			//like E57_JPEG_IMAGE
        int64_t					nImageWidth = 0;
        int64_t					nImageHeight = 0;
        int64_t					nImagesSize = 0;	//Number of bytes in the image
        e57::Image2DType		imageMaskType;		//like E57_PNG_IMAGE_MASK if present
        e57::Image2DType		imageVisualType;	//like E57_JPEG_IMAGE if present

        eReader.GetImage2DSizes(imageIndex, imageProjection, imageType,
            nImageWidth, nImageHeight, nImagesSize, imageMaskType, imageVisualType);

        //Get pixel information off the sphericalRepresentation if imageProjection == E57_SPHERICAL

        int32_t imageHeight = imageHeader.sphericalRepresentation.imageHeight;
        int32_t imageWidth = imageHeader.sphericalRepresentation.imageWidth;
        double pixelHeight = imageHeader.sphericalRepresentation.pixelHeight;
        double pixelWidth = imageHeader.sphericalRepresentation.pixelWidth;

        //Set up buffers

        void* jpegBuffer = new char[nImagesSize];

        //Read the picture data

        eReader.ReadImage2DData(imageIndex, imageProjection, imageType, jpegBuffer, 0, nImageWidth);

        // ... access the picture and decode ...

        //Close and clean up

        delete jpegBuffer;
    }

    write_station_to_bin(station_name, nColumn, nRow, matrix, points, colors, intensity);

    ///////////////////////////////////////////////////////////////////////
    // CATCH THE ERRORS

    //std::cout << " done." << std::endl;
    //std::cout << "Wrinting .bin files...";

    //std::cout << " done." << std::endl;
}

template <class Matrix, class Points, class Colors>
void write_station_to_bin(const std::string &station_name,
                          int nColumn, int nRow,
                          const Matrix &matrix,
                          const std::vector<Points> &points,
                          const std::vector<Colors> &colors,
                          const std::vector<double> &intensity)
{
    enum file_type
    {
        PTS = 0,
        PTX,
        UNSTRUCTURED,
        DOUBLE,
        NONE
    };

    file_type type = file_type::PTX;
    uint64_t nbP = nRow * nColumn;
    uint64_t w = nColumn;
    uint64_t h = nRow;

    std::string matrix_filename = station_name + "_points.tform";
    std::ofstream matrix_file(matrix_filename);
    matrix_file << matrix.m00 << " " << matrix.m01 << " " << matrix.m02 << " " << matrix.m03 << std::endl;
    matrix_file << matrix.m10 << " " << matrix.m11 << " " << matrix.m12 << " " << matrix.m13 << std::endl;
    matrix_file << matrix.m20 << " " << matrix.m21 << " " << matrix.m22 << " " << matrix.m23 << std::endl;
    matrix_file << matrix.m30 << " " << matrix.m31 << " " << matrix.m32 << " " << matrix.m33 << std::endl;
    matrix_file.close();

    std::string out_points = station_name + "_points.bin";
    std::ofstream file_out_points(out_points, std::ios::binary);
    file_out_points.write(reinterpret_cast<char*>(&type), sizeof(file_type));
    file_out_points.write(reinterpret_cast<char*>(&nbP), sizeof(uint64_t));
    file_out_points.write(reinterpret_cast<char*>(&w), sizeof(uint64_t));
    file_out_points.write(reinterpret_cast<char*>(&h), sizeof(uint64_t));
    file_out_points.write(reinterpret_cast<const char*>(points.data()), nbP * sizeof(Points));
    file_out_points.close();

    std::string out_colors = station_name + "_colors.bin";
    std::ofstream file_out_colors(out_colors, std::ios::binary);
    // write first line in colors
    file_out_colors.write(reinterpret_cast<char*>(&nbP), sizeof(uint64_t));
    file_out_colors.write(reinterpret_cast<char*>(&w), sizeof(uint64_t));
    file_out_colors.write(reinterpret_cast<char*>(&h), sizeof(uint64_t));
    file_out_colors.write(reinterpret_cast<const char*>(colors.data()), nbP * sizeof(Colors));
    file_out_colors.close();
}

template <class Matrix, class Points, class Colors>
void write_station_to_ptx(const std::string &station_name,
                          int nColumn, int nRow,
                          const Matrix &matrix,
                          const std::vector<Points> &points,
                          const std::vector<Colors> &colors,
                          const std::vector<double> &intensity)
{
    std::string ptx_fn = station_name + ".ptx";
    std::ofstream ptx_file;
    ptx_file.open(ptx_fn.c_str());
    if (!ptx_file.is_open())
    {
        throw "cloud not open ptx file for writting";
    }
    ptx_file << nColumn << std::endl;
    ptx_file << nRow << std::endl;
    ptx_file << matrix.m03 << " " << matrix.m13 << " " << matrix.m23 << std::endl;
    ptx_file << matrix.m00 << " " << matrix.m01 << " " << matrix.m02 << std::endl;
    ptx_file << matrix.m10 << " " << matrix.m11 << " " << matrix.m12 << std::endl;
    ptx_file << matrix.m20 << " " << matrix.m21 << " " << matrix.m22 << std::endl;
    ptx_file << matrix.m00 << " " << matrix.m01 << " " << matrix.m02 << " " << "0.00000" << std::endl;
    ptx_file << matrix.m10 << " " << matrix.m11 << " " << matrix.m12 << " " << "0.00000" << std::endl;
    ptx_file << matrix.m20 << " " << matrix.m21 << " " << matrix.m22 << " " << "0.00000" << std::endl;
    ptx_file << matrix.m03 << " " << matrix.m13 << " " << matrix.m23 << " " << "1.00000" << std::endl;

    auto itP = points.cbegin();
    auto itI = intensity.cbegin();
    auto itC = colors.cbegin();
    for (size_t idx = 0; idx < points.size(); ++idx, ++itP, ++itI, ++itC)
    {
        ptx_file << itP->x << " " << itP->y << " " << itP->z << " "
                 << *itI << " "
                 << int(itC->r) << " " << int(itC->g) << " " << int(itC->b) << std::endl;
    }

    ptx_file.close();
}

void printNames(int nbParam, char *param[])
{
    // Create a ReaderImpl
    e57::ustring bsFile(param[1]);			//converts Unicode to UTF-8
    e57::Reader	eReader(bsFile);

    ///////////////////////////////////////////////////////////////
    // ACCESSING ROOT DATA

    // Read the root

    e57::E57Root	rootHeader;
    eReader.GetE57Root(rootHeader);

    ///////////////////////////////////////////////////////////////
    // ACCESSING SCAN DATA3D

    //Get the number of scan images available

    int data3DCount = eReader.GetData3DCount();

    //Read the stations' header.
    e57::Data3D		scanHeader;
    for (int scanIndex = 0; scanIndex < data3DCount; ++scanIndex) {
        eReader.ReadData3D(scanIndex, scanHeader);
        std::cout << "scanIndex: " << scanIndex << "; name: " << scanHeader.name
                  << "; description: " << scanHeader.description << std::endl;
    }
}

void getGroupsOBJ(int nbParam, char *param[])
{
    std::ifstream fin(param[1]);
    const int n = 1000;
    char line[n];
    std::map<std::string, std::string> groups;
    while (fin.good())
    {
        fin.getline(line, n);
        if (line[0] == 'g' && line[1] == ' ')
        {
            std::string line_s(line);
            std::string gname = line_s.substr(2);

            fin.getline(line, n);
            line_s = line;
            if (line_s.find("usemtl") == std::string::npos)
            {
                std::cout << "problem" << std::endl;
            }

            std::string matname = line_s.substr(7);
            groups[gname] = matname;
        }
    }
}
