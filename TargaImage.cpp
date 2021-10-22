///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "Globals.h"
#include "TargaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

// constants
const int           RED = 0;                // red channel
const int           GREEN = 1;                // green channel
const int           BLUE = 2;                // blue channel
const unsigned char BACKGROUND[3] = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1; i <= s; i++)
        res = (n - i + 1) * res / i;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
    data = new unsigned char[width * height * 4];
    ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char* d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
        data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image)
{
    width = image.width;
    height = image.height;
    data = NULL;
    if (image.data != NULL) {
        data = new unsigned char[width * height * 4];
        memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
    }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char* rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (!data)
        return NULL;

    // Divide out the alpha
    for (i = 0; i < height; i++)
    {
        int in_offset = i * width * 4;
        int out_offset = i * width * 3;

        for (j = 0; j < width; j++)
        {
            RGBA_To_RGB(data + (in_offset + j * 4), rgb + (out_offset + j * 3));
        }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char* filename)
{
    TargaImage* out_image = Reverse_Rows();

    if (!out_image)
        return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
        cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
        return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char* filename)
{
    unsigned char* temp_data;
    TargaImage* temp_image;
    TargaImage* result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
        width = height = 0;
        return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    unsigned char* RGB_data;
    double* gray_data = new double[width * height];

    RGB_data = To_RGB();

    for (int i = 0, j = 0; i < width * height; ++i, j += 3)
        gray_data[i] = RGB_data[j] * 0.299 + RGB_data[j + 1] * 0.587 + RGB_data[j + 2] * 0.114;

    for (int i = 0, j = 0; i < width * height * 4; i += 4, ++j)
        data[i] = data[i + 1] = data[i + 2] = (double)gray_data[j] / 255 * data[i + 3];

    delete[] RGB_data;
    delete[] gray_data;/**/

    //ClearToBlack();
    return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    unsigned char* RGB_data;
    RGB_data = To_RGB();

    for (int i = 0, j = 0; i < width * height * 4; i += 4, j += 3)
    {
        RGB_data[j] >>= 5;
        RGB_data[j] <<= 5;
        data[i] = (double)RGB_data[j] / 255 * data[i + 3];

        RGB_data[j + 1] >>= 5;
        RGB_data[j + 1] <<= 5;
        data[i + 1] = (double)RGB_data[j + 1] / 255 * data[i + 3];

        RGB_data[j + 2] >>= 6;
        RGB_data[j + 2] <<= 6;
        data[i + 2] = (double)RGB_data[j + 2] / 255 * data[i + 3];
    }

    delete[] RGB_data;

    //ClearToBlack();
    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
    struct Color
    {
        unsigned char r, g, b;
        int num, count = 0;
        static bool compare(Color& i, Color& j)
        {
            return i.count > j.count;
        }
        static bool sequence(Color& i, Color& j)
        {
            return i.num < j.num;
        }
    };

    unsigned char* RGB_data;

    RGB_data = To_RGB();

    for (int i = 0; i < width * height * 3; ++i)
        RGB_data[i] >>= 3;

    Color color[32768];

    for (int i = 0; i < 32768; ++i)
    {
        color[i].num = i;
        color[i].r = i % 32;
        color[i].g = (i % 1024) / 32;
        color[i].b = i / 1024;
    }

    for (int i = 0; i < width * height * 3; i += 3)
    {
        int number;
        number = RGB_data[i] + (RGB_data[i + 1] << 5) + (RGB_data[i + 2] << 10);
        color[number].count++;
    }

    sort(color, color + 32768, Color::compare);

    for (int i = 256; i < 32768; ++i)
    {
        double minDistance = 10000.0;
        double distance = 0.0;
        int index = i;

        for (int j = 0; j < 256; ++j)
        {
            distance = sqrt(pow((double)color[i].r - (double)color[j].r, 2) +
                pow((double)color[i].g - (double)color[j].g, 2) +
                pow((double)color[i].b - (double)color[j].b, 2));

            if (distance < minDistance)
            {
                minDistance = distance;
                index = j;
            }
        }

        color[i].r = color[index].r;
        color[i].g = color[index].g;
        color[i].b = color[index].b;
    }

    sort(color, color + 32768, Color::sequence);

    for (int i = 0, k = 0; i < width * height * 4; i += 4, k += 3)
    {
        int index = RGB_data[k] + (RGB_data[k + 1] << 5) + (RGB_data[k + 2] << 10);

        data[i] = (double)(color[index].r << 3) / 255.0 * data[i + 3];
        data[i + 1] = (double)(color[index].g << 3) / 255.0 * data[i + 3];
        data[i + 2] = (double)(color[index].b << 3) / 255.0 * data[i + 3];
    }

    delete[] RGB_data;

    //ClearToBlack();
    return true;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    To_Grayscale();

    double gray_data = 0.0;

    for (int i = 0; i < width * height * 4; i += 4)
    {
        gray_data = (double)data[i] / (double)data[i + 3];

        if (gray_data < 0.5)
            data[i] = data[i + 1] = data[i + 2] = (double)0 * data[i + 3];
        else
            data[i] = data[i + 1] = data[i + 2] = (double)1 * data[i + 3];
    }

    //ClearToBlack();
    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    double gray_data = 0.0;

    To_Grayscale();

    for (int i = 0; i < width * height * 4; i += 4)
    {
        gray_data = ((double)data[i] / (double)data[i + 3]) + 0.4 * rand() / RAND_MAX - 0.2;
        if (gray_data < 0.5)
            data[i] = data[i + 1] = data[i + 2] = (double)0 * data[i + 3];
        else
            data[i] = data[i + 1] = data[i + 2] = (double)1 * data[i + 3];
    }

    //ClearToBlack();
    return true;
}// Dither_Random


///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
    To_Grayscale();

    int now = 0, next = 0;
    double* gray_data = new double[width * height]();

    for (int i = 0, j = 0; i < width * height * 4; i += 4, ++j)
        gray_data[j] = (double)data[i] / (double)data[i + 3];

    double error = 0.0;

    for (int i = 0; i < height; ++i)
    {
        if (!(i % 2))
        {
            for (int j = 0; j < width; ++j)
            {
                now = j + i * width;
                next = now + width;

                if (gray_data[now] < 0.5)
                {
                    error = gray_data[now];
                    gray_data[now] = 0;
                }
                else
                {
                    error = gray_data[now] - 1.0;
                    gray_data[now] = 1.0;
                }

                if (j == 0)
                {
                    if (i == height - 1)
                        gray_data[now + 1] += error * 0.4375;
                    else {
                        gray_data[now + 1] += error * 0.4375;
                        gray_data[next] += error * 0.3125;
                        gray_data[next + 1] += error * 0.0625;
                    }
                }
                else if (j == width - 1)
                {
                    if (i != height - 1)
                    {
                        gray_data[next] += error * 0.3125;
                        gray_data[next - 1] += error * 0.1875;
                    }
                }
                else
                {
                    if (i == height - 1)
                        gray_data[now + 1] += error * 0.4375;
                    else
                    {
                        gray_data[now + 1] += error * 0.4375;
                        gray_data[next - 1] += error * 0.1875;
                        gray_data[next] += error * 0.3125;
                        gray_data[next + 1] += error * 0.0625;
                    }
                }
            }
        }
        else
        {
            for (int j = width - 1; j >= 0; --j)
            {
                now = j + i * width;
                next = now + width;

                if (gray_data[now] < 0.5)
                {
                    error = gray_data[now];
                    gray_data[now] = 0;
                }
                else
                {
                    error = gray_data[now] - 1.0;
                    gray_data[now] = 1.0;
                }

                if (j == 0)
                {
                    if (i != height - 1)
                    {
                        gray_data[next] += error * 0.3125;
                        gray_data[next + 1] += error * 0.1875;
                    }
                }
                else if (j == width - 1)
                {
                    if (i == height - 1)
                        gray_data[now - 1] += error * 0.4375;
                    else
                    {
                        gray_data[now - 1] += error * 0.4375;
                        gray_data[next] += error * 0.3125;
                        gray_data[next - 1] += error * 0.0625;
                    }
                }
                else
                {
                    if (i == height - 1)
                        gray_data[now - 1] += error * 0.4375;
                    else
                    {
                        gray_data[now - 1] += error * 0.4375;
                        gray_data[next + 1] += error * 0.1875;
                        gray_data[next] += error * 0.3125;
                        gray_data[next - 1] += error * 0.0625;

                    }
                }
            }
        }
    }

    for (int i = 0, j = 0; i < width * height * 4; i += 4, ++j)
        data[i] = data[i + 1] = data[i + 2] = (double)gray_data[j] * data[i + 3];

    //ClearToBlack();
    return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    To_Grayscale();

    double average = 0.0, gray_data = 0.0, threshold = 0.0;
    int hisogram[256];

    for (int i = 0; i < 256; ++i)
        hisogram[i] = 0;

    for (int i = 0; i < width * height * 4; i += 4)
    {
        average += (double)data[i];
        hisogram[data[i]]++;
    }

    average /= (double)(width * height);

    threshold = (1 - average / 255.0) * width * height;

    int k = 0;
    for (; k < 256; ++k)
    {
        threshold -= hisogram[k];
        if (threshold <= 0)
            break;
    }

    if (k == 256)
        threshold = 1.0;
    else
        threshold = (double)k / 255.0;

    for (int i = 0; i < width * height * 4; i += 4)
    {
        gray_data = (double)data[i] / (double)data[i + 3];
        if (gray_data < threshold)
            data[i] = data[i + 1] = data[i + 2] = (double)0 * data[i + 3];
        else
            data[i] = data[i + 1] = data[i + 2] = (double)1 * data[i + 3];
    }

    //ClearToBlack();
    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    double mask[][4] = { {0.7059, 0.3529, 0.5882, 0.2353},
                          {0.0588, 0.9412, 0.8235, 0.4118},
                          {0.4706, 0.7647, 0.8824, 0.1176},
                          {0.1765, 0.5294, 0.2941, 0.6471} };

    int x, y;

    To_Grayscale();

    double gray_data = 0.0;

    for (int i = 0; i < width * height * 4; i += 4)
    {
        x = i / 4 % width;
        y = i / 4 / width;
        gray_data = (double)data[i] / (double)data[i + 3];
        if (gray_data < mask[x % 4][y % 4])
            data[i] = data[i + 1] = data[i + 2] = (double)0 * data[i + 3];
        else
            data[i] = data[i + 1] = data[i + 2] = (double)1 * data[i + 3];

    }

    //ClearToBlack();
    return true;
}// Dither_Cluster


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    unsigned char red_green[8] = { 0, 36, 73, 109, 146, 182, 219, 255 };
    unsigned char blue[4] = { 0, 85, 170, 255 };
    int now = 0, next = 0, index;
    int now_left, now_right, next_left, next_right;
    double error_r = 0, error_g = 0, error_b = 0;
    double distance = 0.0, min_distance = 10000.0;
    unsigned char* temp_data;
    double* RGB_data = new double[width * height * 3]();

    struct Color
    {
        unsigned char r, g, b;
    };

    temp_data = To_RGB();

    for (int i = 0; i < width * height * 3; ++i)
        RGB_data[i] = temp_data[i];

    Color color[256];

    int x = 0;
    for (int i = 0; i < 8; ++i)
        for (int j = 0; j < 8; ++j)
            for (int k = 0; k < 4; ++k)
            {
                color[x].r = red_green[i];
                color[x].g = red_green[j];
                color[x].b = blue[k];
                ++x;
            }

    for (int i = 0; i < height; ++i)
    {
        //row = 0, 2, 4, 6, 8, ...
        if (i % 2 == 0)
        {
            //from 0 to width - 1
            for (int j = 0; j < width; ++j)
            {
                now = (i * width + j) * 3;
                next = now + width * 3;

                now_left = now - 3;
                now_right = now + 3;
                next_left = next - 3;
                next_right = next + 3;

                min_distance = 10000.0;
                distance = 0.0;
                index = 0;
                for (int k = 0; k < 256; ++k)
                {
                    distance = sqrt(pow((double)color[k].r - (double)RGB_data[now], 2) +
                        pow((double)color[k].g - (double)RGB_data[now + 1], 2) +
                        pow((double)color[k].b - (double)RGB_data[now + 2], 2));

                    if (distance < min_distance)
                    {
                        min_distance = distance;
                        index = k;
                    }
                }

                error_r = RGB_data[now] - color[index].r;
                error_g = RGB_data[now + 1] - color[index].g;
                error_b = RGB_data[now + 2] - color[index].b;

                RGB_data[now] = color[index].r;
                RGB_data[now + 1] = color[index].g;
                RGB_data[now + 2] = color[index].b;

                //first pixel of the row
                if (j == 0)
                {
                    if (i == height - 1)
                    {
                        RGB_data[now_right] += error_r * 0.4375;
                        if (RGB_data[now_right] > 255)
                            RGB_data[now_right] = 255;
                        RGB_data[now_right + 1] += error_g * 0.4375;
                        if (RGB_data[now_right + 1] > 255)
                            RGB_data[now_right + 1] = 255;
                        RGB_data[now_right + 2] += error_b * 0.4375;
                        if (RGB_data[now_right + 2] > 255)
                            RGB_data[now_right + 2] = 255;
                    }
                    else
                    {
                        RGB_data[now_right] += error_r * 0.4375;
                        if (RGB_data[now_right] > 255)
                            RGB_data[now_right] = 255;
                        RGB_data[now_right + 1] += error_g * 0.4375;
                        if (RGB_data[now_right + 1] > 255)
                            RGB_data[now_right + 1] = 255;
                        RGB_data[now_right + 2] += error_b * 0.4375;
                        if (RGB_data[now_right + 2] > 255)
                            RGB_data[now_right + 2] = 255;

                        RGB_data[next] += error_r * 0.3125;
                        if (RGB_data[next] > 255)
                            RGB_data[next] = 255;
                        RGB_data[next + 1] += error_g * 0.3125;
                        if (RGB_data[next + 1] > 255)
                            RGB_data[next + 1] = 255;
                        RGB_data[next + 2] += error_b * 0.3125;
                        if (RGB_data[next + 2] > 255)
                            RGB_data[next + 2] = 255;

                        RGB_data[next_right] += error_r * 0.0625;
                        if (RGB_data[next_right] > 255)
                            RGB_data[next_right] = 255;
                        RGB_data[next_right + 1] += error_g * 0.0625;
                        if (RGB_data[next_right + 1] > 255)
                            RGB_data[next_right + 1] = 255;
                        RGB_data[next_right + 2] += error_b * 0.0625;
                        if (RGB_data[next_right + 2] > 255)
                            RGB_data[next_right + 2] = 255;
                    }
                }
                //last pixel of the row
                else if (j == width - 1)
                {
                    if (i != height - 1)
                    {
                        RGB_data[next] += error_r * 0.3125;
                        if (RGB_data[next] > 255)
                            RGB_data[next] = 255;
                        RGB_data[next + 1] += error_g * 0.3125;
                        RGB_data[next + 2] += error_b * 0.3125;

                        RGB_data[next_left] += error_r * 0.1875;
                        RGB_data[next_left + 1] += error_g * 0.1825;
                        RGB_data[next_left + 2] += error_b * 0.1825;
                    }
                }
                //middle pixels of the row
                else
                {
                    if (i == height - 1)
                    {
                        RGB_data[now_right] += error_r * 0.4375;
                        RGB_data[now_right + 1] += error_g * 0.4375;
                        RGB_data[now_right + 2] += error_b * 0.4375;
                    }
                    else
                    {
                        RGB_data[now_right] += error_r * 0.4375;
                        RGB_data[now_right + 1] += error_g * 0.4375;
                        RGB_data[now_right + 2] += error_b * 0.4375;

                        RGB_data[next] += error_r * 0.3125;
                        RGB_data[next + 1] += error_g * 0.3125;
                        RGB_data[next + 2] += error_b * 0.3125;

                        RGB_data[next_right] += error_r * 0.0625;
                        RGB_data[next_right + 1] += error_g * 0.0625;
                        RGB_data[next_right + 2] += error_b * 0.0625;

                        RGB_data[next_left] += error_r * 0.1875;
                        RGB_data[next_left + 1] += error_g * 0.1825;
                        RGB_data[next_left + 2] += error_b * 0.1825;
                    }
                }/**/
            }
        }
        //row = 1, 3, 5, 7, 9, ...
        else
        {
            //from width - 1 to 0
            for (int j = width - 1; j >= 0; --j)
            {
                now = (i * width + j) * 3;
                next = now + width * 3;

                now_left = now - 3;
                now_right = now + 3;
                next_left = next - 3;
                next_right = next + 3;

                min_distance = 10000.0;
                distance = 0.0;
                index = 0;
                for (int k = 0; k < 256; ++k)
                {
                    distance = sqrt(pow((double)color[k].r - (double)RGB_data[now], 2) +
                        pow((double)color[k].g - (double)RGB_data[now + 1], 2) +
                        pow((double)color[k].b - (double)RGB_data[now + 2], 2));

                    if (distance < min_distance)
                    {
                        min_distance = distance;
                        index = k;
                    }
                }

                error_r = RGB_data[now] - color[index].r;
                error_g = RGB_data[now + 1] - color[index].g;
                error_b = RGB_data[now + 2] - color[index].b;

                RGB_data[now] = color[index].r;
                RGB_data[now + 1] = color[index].g;
                RGB_data[now + 2] = color[index].b;

                //last pixel of the row
                if (j == 0)
                {
                    if (i != height - 1)
                    {
                        RGB_data[next] += error_r * 0.3125;
                        RGB_data[next + 1] += error_g * 0.3125;
                        RGB_data[next + 2] += error_b * 0.3125;

                        RGB_data[next_right] += error_r * 0.1875;
                        RGB_data[next_right + 1] += error_g * 0.1875;
                        RGB_data[next_right + 2] += error_b * 0.1875;
                    }
                }
                //first pixel of the row
                else if (j == width - 1)
                {
                    if (i == height - 1)
                    {
                        RGB_data[now_left] += error_r * 0.4375;
                        RGB_data[now_left + 1] += error_g * 0.4375;
                        RGB_data[now_left + 2] += error_b * 0.4375;
                    }
                    else
                    {
                        RGB_data[now_left] += error_r * 0.4375;
                        RGB_data[now_left + 1] += error_g * 0.4375;
                        RGB_data[now_left + 2] += error_b * 0.4375;

                        RGB_data[next] += error_r * 0.3125;
                        RGB_data[next + 1] += error_g * 0.3125;
                        RGB_data[next + 2] += error_b * 0.3125;

                        RGB_data[next_left] += error_r * 0.0625;
                        RGB_data[next_left + 1] += error_g * 0.0625;
                        RGB_data[next_left + 2] += error_b * 0.0625;
                    }
                }
                //middle pixels of the row
                else
                {
                    if (i == height - 1)
                    {
                        RGB_data[now_left] += error_r * 0.4375;
                        RGB_data[now_left + 1] += error_g * 0.4375;
                        RGB_data[now_left + 2] += error_b * 0.4375;
                    }
                    else
                    {
                        RGB_data[now_left] += error_r * 0.4375;
                        RGB_data[now_left + 1] += error_g * 0.4375;
                        RGB_data[now_left + 2] += error_b * 0.4375;

                        RGB_data[next] += error_r * 0.3125;
                        RGB_data[next + 1] += error_g * 0.3125;
                        RGB_data[next + 2] += error_b * 0.3125;

                        RGB_data[next_left] += error_r * 0.0625;
                        RGB_data[next_left + 1] += error_g * 0.0625;
                        RGB_data[next_left + 2] += error_b * 0.0625;

                        RGB_data[next_right] += error_r * 0.1875;
                        RGB_data[next_right + 1] += error_g * 0.1875;
                        RGB_data[next_right + 2] += error_b * 0.1875;
                    }
                }
            }
        }
    }

    for (int i = 0, j = 0; i < width * height * 4; i += 4, j += 3)
    {
        data[i] = (double)RGB_data[j] / 255.0 * data[i + 3];
        data[i + 1] = (double)RGB_data[j + 1] / 255.0 * data[i + 3];
        data[i + 2] = (double)RGB_data[j + 2] / 255.0 * data[i + 3];
    }

    //ClearToBlack();
    return false;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Over: Images not the same size\n";
        return false;
    }

    //ClearToBlack();
    return false;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0; i < width * height * 4; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i + 1] = abs(rgb1[1] - rgb2[1]);
        data[i + 2] = abs(rgb1[2] - rgb2[2]);
        data[i + 3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
    unsigned char* RGB_data;
    unsigned int* temp = new unsigned int[width * height * 3]();

    RGB_data = To_RGB();

    int x, y, pos_x, pos_y;

    for (int i = 0; i < width * height; ++i)
    {
        int now = 0;

        x = i % width;
        y = i / width;

        now = (y * width + x) * 3;

        for (int j = 0; j < 5; ++j)
        {
            pos_x = abs(x - 2 + j);
            if (pos_x > width - 1)
                pos_x = 2 * (width - 1) - pos_x;

            for (int k = 0; k < 5; ++k)
            {
                pos_y = abs(y - 2 + k);
                if (pos_y > height - 1)
                    pos_y = 2 * (height - 1) - pos_y;

                int pos_now = 0;

                pos_now = (pos_y * width + pos_x) * 3;

                temp[now] += RGB_data[pos_now];
                temp[now + 1] += RGB_data[pos_now + 1];
                temp[now + 2] += RGB_data[pos_now + 2];
            }
        }
    }

    for (int i = 0, j = 0; i < width * height * 4; i += 4, j += 3)
    {
        data[i] = ((double)temp[j] * 0.04) / 255.0 * data[i + 3];
        data[i + 1] = ((double)temp[j + 1] * 0.04) / 255.0 * data[i + 3];
        data[i + 2] = ((double)temp[j + 2] * 0.04) / 255.0 * data[i + 3];
    }

    //ClearToBlack();
    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    double mask[5][5] = { { 0.0123, 0.0247, 0.0370, 0.0247, 0.0123 },
                          { 0.0247, 0.0493, 0.0740, 0.0493, 0.0247 },
                          { 0.0370, 0.0740, 0.1111, 0.0740, 0.0370 },
                          { 0.0247, 0.0493, 0.0740, 0.0493, 0.0247 },
                          { 0.0123, 0.0247, 0.0370, 0.0247, 0.0123 } };

    unsigned char* RGB_data;
    double* temp = new double[width * height * 3]();

    RGB_data = To_RGB();

    int x, y, pos_x, pos_y;

    for (int i = 0; i < width * height; ++i)
    {
        int now = 0;

        x = i % width;
        y = i / width;

        now = (y * width + x) * 3;

        for (int j = 0; j < 5; ++j)
        {
            pos_x = abs(x - 2 + j);
            if (pos_x > width - 1)
                pos_x = 2 * (width - 1) - pos_x;

            for (int k = 0; k < 5; ++k)
            {
                pos_y = abs(y - 2 + k);
                if (pos_y > height - 1)
                    pos_y = 2 * (height - 1) - pos_y;

                int pos_now = 0;

                pos_now = (pos_y * width + pos_x) * 3;

                temp[now] += RGB_data[pos_now] * mask[j][k];
                temp[now + 1] += RGB_data[pos_now + 1] * mask[j][k];
                temp[now + 2] += RGB_data[pos_now + 2] * mask[j][k];
            }
        }
    }

    for (int i = 0, j = 0; i < width * height * 4; i += 4, j += 3)
    {
        data[i] = (double)temp[j] / 255.0 * data[i + 3];
        data[i + 1] = (double)temp[j + 1] / 255.0 * data[i + 3];
        data[i + 2] = (double)temp[j + 2] / 255.0 * data[i + 3];
    }
    //ClearToBlack();
    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    double mask[5][5] = { { 0.0039, 0.0156, 0.0234, 0.0156, 0.0039 },
                          { 0.0156, 0.0625, 0.0938, 0.0625, 0.0156 },
                          { 0.0234, 0.0938, 0.1406, 0.0938, 0.0234 },
                          { 0.0156, 0.0625, 0.0938, 0.0625, 0.0156 },
                          { 0.0039, 0.0156, 0.0234, 0.0156, 0.0039 } };

    unsigned char* RGB_data;
    double* temp = new double[width * height * 3]();

    RGB_data = To_RGB();

    int x, y, pos_x, pos_y;

    for (int i = 0; i < width * height; ++i)
    {
        int now = 0;

        x = i % width;
        y = i / width;

        now = (y * width + x) * 3;

        for (int j = 0; j < 5; ++j)
        {
            pos_x = abs(x - 2 + j);
            if (pos_x > width - 1)
                pos_x = 2 * (width - 1) - pos_x;

            for (int k = 0; k < 5; ++k)
            {
                pos_y = abs(y - 2 + k);
                if (pos_y > height - 1)
                    pos_y = 2 * (height - 1) - pos_y;

                int pos_now = 0;

                pos_now = (pos_y * width + pos_x) * 3;

                temp[now] += RGB_data[pos_now] * mask[j][k];
                temp[now + 1] += RGB_data[pos_now + 1] * mask[j][k];
                temp[now + 2] += RGB_data[pos_now + 2] * mask[j][k];
            }
        }
    }

    for (int i = 0, j = 0; i < width * height * 4; i += 4, j += 3)
    {
        data[i] = (double)temp[j] / 255.0 * data[i + 3];
        data[i + 1] = (double)temp[j + 1] / 255.0 * data[i + 3];
        data[i + 2] = (double)temp[j + 2] / 255.0 * data[i + 3];
    }

    //ClearToBlack();
    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N(unsigned int N)
{
    int n = 1;
    int* pascal = new int[N]();
    pascal[0] = 1;

    for (int i = 1; i < N; ++i)
        pascal[i] = pascal[i - 1] * (N - i + 1) / i;

    pascal[N - 1] = 1;

    double** mask;

    mask = new double* [N]();

    for (int i = 0; i < N; ++i)
        mask[i] = new double[N]();

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            mask[i][j] = (double) pascal[i] * pascal[j] / pow(N, 2);

    unsigned char* RGB_data;
    double* temp = new double[width * height * 3]();

    RGB_data = To_RGB();

    int x, y, pos_x, pos_y;

    for (int i = 0; i < width * height; ++i)
    {
        int now = 0;

        x = i % width;
        y = i / width;

        now = (y * width + x) * 3;

        for (int j = 0; j < 5; ++j)
        {
            pos_x = abs(x - 2 + j);
            if (pos_x > width - 1)
                pos_x = 2 * (width - 1) - pos_x;

            for (int k = 0; k < 5; ++k)
            {
                pos_y = abs(y - 2 + k);
                if (pos_y > height - 1)
                    pos_y = 2 * (height - 1) - pos_y;

                int pos_now = 0;

                pos_now = (pos_y * width + pos_x) * 3;

                temp[now] += RGB_data[pos_now] * mask[j][k];
                temp[now + 1] += RGB_data[pos_now + 1] * mask[j][k];
                temp[now + 2] += RGB_data[pos_now + 2] * mask[j][k];
            }
        }
    }

    for (int i = 0, j = 0; i < width * height * 4; i += 4, j += 3)
    {
        data[i] = (double)temp[j] / 255.0 * data[i + 3];
        data[i + 1] = (double)temp[j + 1] / 255.0 * data[i + 3];
        data[i + 2] = (double)temp[j + 2] / 255.0 * data[i + 3];
    }

    //ClearToBlack();
    return false;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    double mask[5][5] = { { -0.0039, -0.0156, -0.0234, -0.0156, -0.0039 },
                          { -0.0156, -0.0625, -0.0938, -0.0625, -0.0156 },
                          { -0.0234, -0.0938,  0.8594, -0.0938, -0.0234 },
                          { -0.0156, -0.0625, -0.0938, -0.0625, -0.0156 },
                          { -0.0039, -0.0156, -0.0234, -0.0156, -0.0039 } };

    unsigned char* RGB_data;
    double* temp = new double[width * height * 3]();

    RGB_data = To_RGB();

    int x, y, pos_x, pos_y;

    for (int i = 0; i < width * height; ++i)
    {
        int now = 0;

        x = i % width;
        y = i / width;

        now = (y * width + x) * 3;

        for (int j = 0; j < 5; ++j)
        {
            pos_x = abs(x - 2 + j);
            if (pos_x > width - 1)
                pos_x = 2 * (width - 1) - pos_x;

            for (int k = 0; k < 5; ++k)
            {
                pos_y = abs(y - 2 + k);
                if (pos_y > height - 1)
                    pos_y = 2 * (height - 1) - pos_y;

                int pos_now = 0;

                pos_now = (pos_y * width + pos_x) * 3;

                temp[now] += RGB_data[pos_now] * mask[j][k];
                temp[now + 1] += RGB_data[pos_now + 1] * mask[j][k];
                temp[now + 2] += RGB_data[pos_now + 2] * mask[j][k];
            }
        }
    }

    for (int i = 0, j = 0; i < width * height * 4; i += 4, j += 3)
    {
        if (temp[j] < 0)
            temp[j] = 0;
        if (temp[j + 1] < 0)
            temp[j + 1] = 0;
        if (temp[j + 2] < 0)
            temp[j + 2] = 0;

        data[i] = (double)temp[j] / 255.0 * data[i + 3];
        data[i + 1] = (double)temp[j + 1] / 255.0 * data[i + 3];
        data[i + 2] = (double)temp[j + 2] / 255.0 * data[i + 3];
    }
    //ClearToBlack();
    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    double mask[5][5] = { { -0.0039, -0.0156, -0.0234, -0.0156, -0.0039 },
                          { -0.0156, -0.0625, -0.0938, -0.0625, -0.0156 },
                          { -0.0234, -0.0938,  1.8594, -0.0938, -0.0234 },
                          { -0.0156, -0.0625, -0.0938, -0.0625, -0.0156 },
                          { -0.0039, -0.0156, -0.0234, -0.0156, -0.0039 } };

    unsigned char* RGB_data;
    double* temp = new double[width * height * 3]();

    RGB_data = To_RGB();

    int x, y, pos_x, pos_y;

    for (int i = 0; i < width * height; ++i)
    {
        int now = 0;

        x = i % width;
        y = i / width;

        now = (y * width + x) * 3;

        for (int j = 0; j < 5; ++j)
        {
            pos_x = abs(x - 2 + j);
            if (pos_x > width - 1)
                pos_x = 2 * (width - 1) - pos_x;

            for (int k = 0; k < 5; ++k)
            {
                pos_y = abs(y - 2 + k);
                if (pos_y > height - 1)
                    pos_y = 2 * (height - 1) - pos_y;

                int pos_now = 0;

                pos_now = (pos_y * width + pos_x) * 3;

                temp[now] += RGB_data[pos_now] * mask[j][k];
                temp[now + 1] += RGB_data[pos_now + 1] * mask[j][k];
                temp[now + 2] += RGB_data[pos_now + 2] * mask[j][k];
            }
        }
    }

    for (int i = 0, j = 0; i < width * height * 4; i += 4, j += 3)
    {
        if (temp[j] < 0)
            temp[j] = 0;
        if (temp[j + 1] < 0)
            temp[j + 1] = 0;
        if (temp[j + 2] < 0)
            temp[j + 2] = 0;

        if (temp[j] > 255)
            temp[j] = 255;
        if (temp[j + 1] > 255)
            temp[j + 1] = 255;
        if (temp[j + 2] > 255)
            temp[j + 2] = 255;

        data[i] = (double)temp[j] / 255.0 * data[i + 3];
        data[i + 1] = (double)temp[j + 1] / 255.0 * data[i + 3];
        data[i + 2] = (double)temp[j + 2] / 255.0 * data[i + 3];
    }
    //ClearToBlack();
    return true;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{

    //ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{

    double mask[3][3] = { { 0.0625, 0.1250, 0.0625 },
                          { 0.1250, 0.2500, 0.1250 },
                          { 0.0625, 0.1250, 0.0625} };

    unsigned char* RGB_data;

    int new_width = width / 2, new_height = height / 2;

    double* temp = new double[new_width * new_height * 3]();

    RGB_data = To_RGB();

    int x, y, pos_x, pos_y;

    for (int i = 0; i < new_height; ++i)
    {
        for (int j = 0; j < new_width; ++j)
        {
            int now = 0;

            x = j * 2;
            y = i * 2;

            now = (i * new_width + j) * 3;

            for (int k = 0; k < 3; ++k)
            {
                pos_x = abs(x - 1 + k);
                if (pos_x > width - 1)
                    pos_x = 2 * (width - 1) - pos_x;

                for (int m = 0; m < 3; ++m)
                {
                    pos_y = abs(y - 1 + m);
                    if (pos_y > height - 1)
                        pos_y = 2 * (height - 1) - pos_y;

                    int pos_now = 0;

                    pos_now = (pos_y * width + pos_x) * 3;

                    temp[now] += RGB_data[pos_now] * mask[k][m];
                    temp[now + 1] += RGB_data[pos_now + 1] * mask[k][m];
                    temp[now + 2] += RGB_data[pos_now + 2] * mask[k][m];
                }
            }
        }
    }

    width = new_width;
    height = new_height;

    delete[] data;

    for (int i = 0; i < new_width * new_height * 3; ++i)
        if (temp[i] > 255)
            temp[i] = 255;
        else if (temp[i] < 0)
            temp[i] = 0;

    data = new unsigned char[new_width * new_height * 4]();

    for (int i = 0, j = 0; i < new_width * new_height * 4; i += 4, j += 3)
    {
        data[i] = (double)temp[j];
        data[i + 1] = (double)temp[j + 1];
        data[i + 2] = (double)temp[j + 2];
        data[i + 3] = 255;
    }

    //ClearToBlack();
    return true;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    double maskee[3][3] = { { 0.0625, 0.1250, 0.0625 },
                            { 0.1250, 0.2500, 0.1250 },
                            { 0.0625, 0.1250, 0.0625 } };
    double maskoo[4][4] = { { 0.0156, 0.0469, 0.0469, 0.0156},
                            { 0.0469, 0.1406, 0.1406, 0.0469},
                            { 0.0469, 0.1406, 0.1406, 0.0469},
                            { 0.0156, 0.0469, 0.0469, 0.0156} };
    double maskoe[4][3] = { { 0.0313, 0.0625, 0.0313 },
                            { 0.0938, 0.1875, 0.0938 },
                            { 0.0938, 0.1875, 0.0938 },
                            { 0.0313, 0.0625, 0.0313} };
    double maskeo[3][4] = { { 0.0313, 0.0938, 0.0938, 0.0313 },
                            { 0.0625, 0.1875, 0.1875, 0.0938 },
                            { 0.0313, 0.0938, 0.0938, 0.0313 } };

    unsigned char* RGB_data;

    int new_width = width * 2, new_height = height * 2;

    double* temp = new double[new_width * new_height * 3]();

    RGB_data = To_RGB();

    double x, y;
    int pos_x, pos_y;

    for (int j = 0; j < new_height; ++j)
    {
        for (int i = 0; i < new_width; ++i)
        {
            int now = 0;

            x = (double) i / 2;
            y = (double) j / 2;

            now = (j * new_width + i) * 3;

            if (!(i % 2) && !(j % 2))
            {
                for (int m = 0; m < 3; ++m)
                {
                    pos_x = abs(x - 1 + m);
                    if (pos_x > width - 1)
                        pos_x = 2 * (width - 1) - pos_x;
                    for (int n = 0; n < 3; ++n)
                    {
                        pos_y = abs(y - 1 + n);
                        if (pos_y > height - 1)
                            pos_y = 2 * (height - 1) - pos_y;

                        int pos_now = 0;

                        pos_now = (pos_y * width + pos_x) * 3;

                        temp[now] += RGB_data[pos_now] * maskee[m][n];
                        temp[now + 1] += RGB_data[pos_now + 1] * maskee[m][n];
                        temp[now + 2] += RGB_data[pos_now + 2] * maskee[m][n];
;                   }
                }
            }
            else if (!(i % 2) && (j % 2))
            {
                for (int m = 0; m < 3; ++m)
                {
                    pos_x = abs((double) x - 1 + m);
                    if (pos_x > width - 1)
                        pos_x = 2 * (width - 1) - pos_x;
                    for (int n = 0; n < 2; ++n)
                    {
                        pos_y = abs((double) y - 1 + n);
                        if (pos_y > height - 1)
                            pos_y = 2 * (height - 1) - pos_y;

                        int pos_now = 0;

                        pos_now = (pos_y * width + pos_x) * 3;

                        temp[now] += RGB_data[pos_now] * maskoe[m][n];
                        temp[now + 1] += RGB_data[pos_now + 1] * maskoe[m][n];
                        temp[now + 2] += RGB_data[pos_now + 2] * maskoe[m][n];
                        ;
                    }
                }
            }
            else if ((i % 2) && !(j % 2))
            {
                for (int m = 0; m < 2; ++m)
                {
                    pos_x = abs((double) x - 1 + m);
                    if (pos_x > width - 1)
                        pos_x = 2 * (width - 1) - pos_x;
                    for (int n = 0; n < 3; ++n)
                    {
                        pos_y = abs(y - 1 + n);
                        if (pos_y > height - 1)
                            pos_y = 2 * (height - 1) - pos_y;

                        int pos_now = 0;

                        pos_now = (pos_y * width + pos_x) * 3;

                        temp[now] += RGB_data[pos_now] * maskeo[m][n];
                        temp[now + 1] += RGB_data[pos_now + 1] * maskeo[m][n];
                        temp[now + 2] += RGB_data[pos_now + 2] * maskeo[m][n];
                    }
                }
            }
            else
            {
                for (int m = 0; m < 3; ++m)
                {
                    pos_x = abs((double) x - 1 + m);
                    if (pos_x > width - 1)
                        pos_x = 2 * (width - 1) - pos_x;
                    for (int n = 0; n < 3; ++n)
                    {
                        pos_y = abs((double) y - 1 + n);
                        if (pos_y > height - 1)
                            pos_y = 2 * (height - 1) - pos_y;

                        int pos_now = 0;

                        pos_now = (pos_y * width + pos_x) * 3;

                        temp[now] += RGB_data[pos_now] * maskoo[m][n];
                        temp[now + 1] += RGB_data[pos_now + 1] * maskoo[m][n];
                        temp[now + 2] += RGB_data[pos_now + 2] * maskoo[m][n];
                    }
                }
            }
        }
    }

    width = new_width;
    height = new_height;

    delete[] data;

    for (int i = 0; i < new_width * new_height * 3; ++i)
        if (temp[i] > 255)
            temp[i] = 255;
        else if (temp[i] < 0)
            temp[i] = 0;

    data = new unsigned char[new_width * new_height * 4]();

    for (int i = 0, j = 0; i < new_width * new_height * 4; i += 4, j += 3)
    {
        data[i] = (double)temp[j];
        data[i + 1] = (double)temp[j + 1];
        data[i + 2] = (double)temp[j + 2];
        data[i + 3] = 255;
    }

    //ClearToBlack();
    return false;
}// Double_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char* rgba, unsigned char* rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {

        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
        float	alpha_scale = (float)255 / (float)alpha;
        int	val;
        int	i;

        for (i = 0; i < 3; i++)
        {
            val = (int)floor(rgba[i] * alpha_scale);
            if (val < 0)
                rgb[i] = 0;
            else if (val > 255)
                rgb[i] = 255;
            else
                rgb[i] = val;
        }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char* dest = new unsigned char[width * height * 4];
    TargaImage* result;
    int 	        i, j;

    if (!data)
        return NULL;

    for (i = 0; i < height; i++)
    {
        int in_offset = (height - i - 1) * width * 4;
        int out_offset = i * width * 4;

        for (j = 0; j < width; j++)
        {
            dest[out_offset + j * 4] = data[in_offset + j * 4];
            dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
            dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
            dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
    int radius_squared = (int)s.radius * (int)s.radius;
    for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
        for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
            int x_loc = (int)s.x + x_off;
            int y_loc = (int)s.y + y_off;
            // are we inside the circle, and inside the image?
            if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
                int dist_squared = x_off * x_off + y_off * y_off;
                if (dist_squared <= radius_squared) {
                    data[(y_loc * width + x_loc) * 4 + 0] = s.r;
                    data[(y_loc * width + x_loc) * 4 + 1] = s.g;
                    data[(y_loc * width + x_loc) * 4 + 2] = s.b;
                    data[(y_loc * width + x_loc) * 4 + 3] = s.a;
                }
                else if (dist_squared == radius_squared + 1) {
                    data[(y_loc * width + x_loc) * 4 + 0] =
                        (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
                    data[(y_loc * width + x_loc) * 4 + 1] =
                        (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
                    data[(y_loc * width + x_loc) * 4 + 2] =
                        (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
                    data[(y_loc * width + x_loc) * 4 + 3] =
                        (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
                }
            }
        }
    }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
    unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
    radius(iradius), x(ix), y(iy), r(ir), g(ig), b(ib), a(ia)
{
}

