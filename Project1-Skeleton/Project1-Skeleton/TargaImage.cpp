///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//                                              Modified:   Feng Liu
//                                              Date:       Winter 2011
//                                              Why:        Change the library file 
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
#include <time.h>
#include <random>
#include <valarray>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

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
TargaImage::TargaImage(int w, int h, unsigned char *d)
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
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
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
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
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
	//To_RGB(); - not doing this because it will destroy alpha
	int i, grey = 0, size = width * height * 4;
	//+4 to leave alpha alone
	for(i = 0; i < size; i += 4)
	{
			grey = data[i + RED]*0.299+data[i + GREEN]*0.587+data[i + BLUE]*0.114;
			data[i + RED] = data[i + GREEN] = data[i + BLUE] = grey;
	}
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
	int i, size = width * height * 4;
	for(i = 0; i < size; i +=4)
	{
		data[i+RED] /= 32;
		data[i+GREEN] /= 32;
		data[i+BLUE] /= 64;
	}
	for(i = 0; i < size; i += 4)
	{
		data[i+RED] *= 32;
		data[i+GREEN] *= 32;
		data[i+BLUE] *= 64;
	}
	return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
struct colorHist{
	colorHist(){ r = 0; g = 0; b = 0; /*color = 0;*/ count = 0; }
	//long color;
	int r ;
	int g;
	int b;
	int count;
};

bool sortingFunc(const colorHist & obj1, const colorHist & obj2){
	return (obj1.count < obj2.count);
}
bool TargaImage::Quant_Populosity()
{
	int i, j = 0, size = width * height * 4, histSize = 32*32*32+1; //total possible colors = 32*32*32m(a little bigger because one was outside range), 256 colors * 4 (rgba) for final image
	long startIndex = 0, big = 5000000; //set proximity to some arbitrarily large number to make sure it is properly initialized
	long colorIndex = 0, colorToMap = 0, toCompare = 0;
	double currProx = 0, proximity = big;
	float scalingFactor = 8.0;
	unsigned char r,g,b;
	vector<colorHist> colors(histSize);

	for(i = 0; i < size; i +=4)
	{
		/*data[i+RED] = (unsigned char)(data[i+RED] / scalingFactor);								//scale each color down to 5 bits (scaling factor is remainder of largest possible, so 255, divided by 32)  
		data[i+GREEN] = (unsigned char)(data[i+GREEN] / scalingFactor);
		data[i+BLUE] = (unsigned char)(data[i+BLUE] / scalingFactor);
		colorIndex = (data[i+RED] + (((unsigned int)data[i+GREEN] << 5) + ((unsigned int)data[i+BLUE] << 10)));*/		//store bits of pixel in long, ignore alpha
		r = (unsigned char)(data[i+RED]/scalingFactor);
		g = (unsigned char)(data[i+GREEN]/scalingFactor);
		b = (unsigned char)(data[i+BLUE]/scalingFactor);
		colorIndex = (r + (((unsigned int)g << 5) + ((unsigned int)b << 10)));
		assert(colorIndex < 32768);
		colors[colorIndex].r = (int)data[i+RED];
		colors[colorIndex].g = (int)data[i+GREEN];
		colors[colorIndex].b = (int)data[i+BLUE];
		++colors[colorIndex].count;
		//colors[colorIndex].color = colorIndex;
		//++colors[colorIndex].count;
	}

	sort(colors.begin(),colors.end(),sortingFunc);										//sort colors by number of times they were found in histogram
	startIndex = histSize - 256;														//256 most popular colors will be checked against exisiting colors in array
	for(i = 0; i < size; i+=4)															//Check each color in original data set against most popular colors in histogram
	{
		proximity = big;																//reset the proximity each time we loop
		//toCompare = (data[i+RED] + (((unsigned int)data[i+GREEN] << 5) + ((unsigned int)data[i+BLUE] << 10)));		//convert current color to bits
		for(j = startIndex; j < histSize; ++j)
		{
			/*currProx = abs(colors[j].color - toCompare);								//determine whether this color in most popular colors array is closest to current color
			if(currProx <= proximity) {
				proximity = currProx;
				colorToMap = colors[j].color;											
			}*/
			currProx = sqrt(((data[i+RED]-colors[j].r)*(data[i+RED]-colors[j].r))+((data[i+GREEN]-colors[j].g)*(data[i+GREEN]-colors[j].g))+((data[i+BLUE]-colors[j].b)*(data[i+BLUE]-colors[j].b)));
			if(currProx < proximity){
				proximity = currProx;
				colorToMap = j;
			}
		}
		/*data[i+RED] = (unsigned char)floor((colorToMap & 0xff)*scalingFactor);					 //convert back to 24 bit unsigned char so it will display properly
		data[i+GREEN] = (unsigned char)floor(((colorToMap >> 5) & 0xff)* scalingFactor); 
		data[i+BLUE] = (unsigned char)floor(((colorToMap >> 10) & 0xff) * scalingFactor);*/
		data[i+RED] = (unsigned char)colors[colorToMap].r;
		data[i+GREEN] = (unsigned char)colors[colorToMap].g;
		data[i+BLUE] = (unsigned char)colors[colorToMap].b;
	}
	/*for(i = 0; i < size; i+=4)
	{
		data[i+RED] *= 32;
		data[i+GREEN] *= 32;
		data[i+BLUE] *= 64;
	}*/
    return true;
}// Quant_Populosity
/*r = (colors[j].color & 0xff);
			g = ((colors[j].color >> 5) & 0xff);
			b = ((colors[j].color >> 10) & 0xff);
			currProx = sqrt(((data[i+RED]-r)*(data[i+RED]-r)) + ((data[i+GREEN]-g)*(data[i+GREEN]-g)) + ((data[i+BLUE]-b)*(data[i+BLUE]-b)));*/

///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
	int i, size = width * height * 4, newColor;
	float color;
	To_Grayscale();
	for(i=0; i<size; i+=4)
	{
		color = data[i]/(float)256;
		if(color < 0.5)
			newColor = 0;
		else
			newColor = 255;
		data[i+RED] = data[i+GREEN] = data[i+BLUE] = newColor;
	}
    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
//use random library instead of rand
	float arr[42] = {0};
	int i, size = width * height * 4, newColor, index = 0;
	float color, start = -0.2;
	for(i = 0; i < 42; ++i){
		arr[i] = start;
		start += 0.01;
	}
	To_Grayscale();
	for(i = 0; i < size; i += 4)
	{
		color = data[i]/(float)256;
		index = rand()%42;
		color += arr[index];
		if(color < 0.5)
			newColor = 0;
		else
			newColor = 255;
		data[i] = data[i+GREEN] = data[i+BLUE] = newColor;
	}
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
	int i = 0, j = 0, realWidth = width * 4, size = width * height * 4, count = 0, toMove = 0;
	//int right = 0, down = 0, downR = 0, downL = 0;
	To_Grayscale();
	float * output = new float [size];
	for(i = 0; i < size; i+= 4)	{				//convert to float in range 0-1, ignore alpha
		output[i] = output[i+GREEN] = output[i+BLUE] = (data[i]/(float)256);
		//output[i+3] = (float)data[i+3];
	}
	for(i = 0; i < height; ++i)
	{
		for(j = 0; j < realWidth ; j+=4)
		{
			if(i%2==0)							//we are at an even row, so we should go from left to right
				leftRight(i,j,output);
			else
				rightLeft(i,j,output);			//odd row, so we go right to left
		}
	}
	for(i = 0; i < height; ++i)
	{
		for(j = 0; j < realWidth; j+=4)
		{
			data[i*realWidth+j] = (unsigned char)output[i*realWidth+j];
			data[i*realWidth+j+GREEN] = (unsigned char)output[i*realWidth+j+GREEN];
			data[i*realWidth+j+BLUE] = (unsigned char)output[i*realWidth+j+BLUE];
			//data[i*realWidth+j+3] = (unsigned char)output[i*realWidth+j+3];
		}
	}
	delete [] output;
    return true;
}// Dither_FS

void TargaImage::leftRight(int i, int j, float * output) 
{ 
 	int right, down, downL, downR, realWidth = width * 4; 
 	right = (i * realWidth) + j + 1;//4; 
 	if(j == realWidth-4)  
 		right = -7;			//set to some arbitrary invalid index, this will fail check in function that calculates error 
 	down = (i+1)* realWidth + j; 
 	downL = (i+1)*realWidth+j-4; 
 	downR = (i+1)*realWidth+j+4; 
 	if(i+1 >= height){ 
 		down = -7; 
 		downL = -7; 
 		downR = -7; 
 	} 
 	if(j+4 == realWidth) 
 		downR = -7; 
 	if(j == 0) 
 		downL = -7; 
 	int currIndex = i * realWidth + j; 
 	calcThresh(currIndex,output,right,down,downL,downR); 
} 
  


void TargaImage::rightLeft(int i, int j, float * output) 
{ 
	int right, down, downL, downR, realWidth = width * 4, size = width * height * 4; 
 	right = (i * realWidth) + (realWidth-j) -1; //-4; 
	if(j ==0)  
		right = -7; 
	down = (i+1)*realWidth+(realWidth-j); 
	downL = (i+1)*realWidth+(realWidth-j)-4; 
	downR = (i+1)*realWidth+(realWidth-j)+4; 
	if(i+1 >= height){ 
		down = -7; 
		downL = -7; 
		downR = -7; 
	} 
	if(j == 0) 
		downL = -7; 
	if(j == realWidth-4)
		downR = -7; 
	int currIndex = i * realWidth + j; 
	calcThresh(currIndex,output,right,down,downL,downR); 
}

void TargaImage::calcThresh(int i, float * output, int over, int below, int downL, int downR) //I think I need some logic in here if data[i,j] < 0 or > 1... if so, "ignore", so return?
{
		double error = 0, down = 5.0/16.0, leftDown = 3.0/16.0, right = 7.0/16.0, rightDown = 1.0/16.0;
		int t = 0, result = 0; //TODO: test all against 0, if >= 0 calculate error, else ignore.
		bool lessThan = output[i] < 0.5;
		if(lessThan) t = 0;
		else
			t = 1;
		error = output[i]-t;
		if(below >= 0){
		output[below]+= (down*error);
		output[below+GREEN] = output[below+BLUE] = output[below+3] = output[below];
		}
		if(downR >= 0){
		output[downR] += (rightDown*error); //downOver is right on even rows, left on odd rows
		output[downR+GREEN] = output[downR+BLUE] = output[downR+3] = output[downR];
		}
		if(over >= 0){
		output[over] += (right * error); 
		output[over+GREEN] = output[over+BLUE] = output[over+3] = output[over];
		}
		if(downL >= 0){
			output[downL] += (leftDown * error);
			output[downL+GREEN] = output[downL+BLUE] = output[downL+3] = output[downL];
		}
		if(lessThan)
			result = output[i] = output[i+GREEN] = output[i+BLUE] = output[i+3] = 0;
		else
			result = output[i] = output[i+GREEN] = output[i+BLUE] = output[i+3] = 255;

		//cout << "Pixel " << i << ": " << result << endl;

}

///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Dither_Bright()
{
	int i = 0, count = 0, size = width * height * 4, numPixels = width * height, newColor;
	float sum = 0, average = 0, thresh = 0;
	To_Grayscale();
	vector<unsigned char> copy;
	copy.reserve(numPixels);	//reserve memory for size numPixels
	for(i = 0; i < size; i+=4) { //copy pixel info from data to copy vector
		copy.push_back(data[i]);
		sum += data[i];
	}
	average = sum/(float)(numPixels);
	sort(copy.begin(), copy.end());
	int index = (1-average/256)*numPixels;
	thresh = copy[index];
	for(i = 0; i <size; i+=4) {
		if(data[i] < thresh)
			newColor = 0;
		else
			newColor = 255;
		data[i] = data[i+GREEN] = data[i+BLUE] = newColor;
	}
	//do i need to free vector memory - if so how
    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    ClearToBlack();
    return false;
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
    ClearToBlack();
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
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

    ClearToBlack();
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

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
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
	float filter[25] = {0};
	int i = 0, size = width * height * 4;
    ClearToBlack();
    return false;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    ClearToBlack();
    return false;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
	int i = 0, j = 0, filterSize = 5, sizeRGBA = width * height * 4, rowOffset = i * width * 4;
	double sum = 0.0;
	unsigned char * newData = new unsigned char[sizeRGBA];

	for(i = 0; i < sizeRGBA; ++i)		//initialize entire array to 0 to start
		newData[i] = 0;

	double * mask = generateGaussFilter(filterSize);

	//apply mask to image
	for (i = 2 ; i < height-2 ; ++i)	//ignore edges
    {
		rowOffset = i * width * 4;
		for (j = 2 ; j < width-2 ; ++j)	//ignore edges
			applyFilter(mask, filterSize, newData, rowOffset, i, j);
    }
	for (i = 2; i < height - 2; ++i)
	{
		for (j = 2; j < width - 2; ++j)
		{
			rowOffset = i * width * 4;
			data[rowOffset + j * 4] = newData[rowOffset + j * 4];
			data[(rowOffset + j * 4) + GREEN] = newData[(rowOffset + j * 4) + GREEN];
			data [(rowOffset + j * 4) + BLUE] = newData[(rowOffset + j * 4) + BLUE];
		}
	}
	delete [] mask;
	delete [] newData;
	return true;
}// Filter_Gaussian

void TargaImage::applyFilter(double * mask, int filterSize, unsigned char * newData, int rowOffset, int i, int j)
{
	int sum = 0, color = 0, halfFilterSize = filterSize/2, indexI = 0, indexJ = 0, u = 0, v = 0;
	while(color < 3){
		sum = 0;
		for(int u = -halfFilterSize; u <= halfFilterSize; ++u){
			for(v = -halfFilterSize; v <= halfFilterSize; ++v){
				indexI = i+u;//ADDED THIS -- edge detection for gaussian filter
				indexJ = j+v;
				if(indexI < 0)
					indexI = 0;
				else if(indexI >= height)
					indexI = height - 1;
				if(indexJ >= width)
					indexJ = width-1;
				else if(indexJ < 0)
					indexJ = 0;
				//sum += data[(i+u)*width*4 + (j+v)*4 + color]*mask[(u+halfFilterSize)*filterSize+v+halfFilterSize];	//apply to all red in scope of filter first, then green, then blue
				sum += data[indexI*width*4 + indexJ*4 + color]*mask[(u+halfFilterSize)*filterSize+v+halfFilterSize];
			}
		}
		newData[(rowOffset + j*4) + color] = sum;
		++color;
	}
}

double* TargaImage::generateGaussFilter(int n)
{
	int size = n*n;
	double * mask = new double[size], first = 0.0, second = 0.0;
	int i = 0, j = 0;
	double sum = 0.0;
	for(i = 0; i < n; ++i)
	{
		for(j = 0; j < n; ++j){
			first = Binomial(n-1,i);
			second = Binomial(n-1,j);
			mask[i*n+j] = first*second;
			sum += first*second;
		}
	}
	for(i = 0; i < size; ++i)
	{
		mask[i]/=sum;
	}
	return mask;

}
///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N( unsigned int N )
{
	//The issue is that the filter does not recognize the edges, and so the image gets darker because it is being summed with black edges (I think)
   	int i = 0, j = 0, filterSize = N, halfFilterWidth = N/2, sizeRGBA = width * height * 4, rowOffset;
	double sum = 0.0;
	unsigned char * newData = new unsigned char[sizeRGBA];
	for(i = 0; i < sizeRGBA; ++i)
		newData[i] = 0;
	double * mask = generateGaussFilter(filterSize);
	//apply mask to image
	for (int i = 0; i < height; ++i)
    {
		rowOffset = i * width * 4;
		for (int j = 0; j < width; ++j)
		{
			applyFilter(mask, filterSize, newData, rowOffset, i, j);
		}
    }
	for (int i = 0; i < height; ++i)
	{
		rowOffset = i * width * 4;
		for (int j = 0; j < width; ++j)
		{
			data[rowOffset + j * 4] = newData[rowOffset + j * 4];
			data[(rowOffset + j * 4) + GREEN] = newData[(rowOffset + j * 4) + GREEN];
			data[(rowOffset + j * 4) + BLUE] = newData[(rowOffset + j * 4) + BLUE];
		}
	}
	delete [] mask;
	delete [] newData;
	return true;
   return true;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    ClearToBlack();
    return false;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    ClearToBlack();
    return false;
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
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    ClearToBlack();
    return false;
}// Half_Size


///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    ClearToBlack();
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
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
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

	    for (i = 0 ; i < 3 ; i++)
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
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
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
            } else if (dist_squared == radius_squared + 1) {
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

void TargaImage::To_RGBA(unsigned char * newData)
{
	int i = 0, j = 0;
	for(i = 0; i < height; ++i)
	{
		int out_offset = i * width * 4;
		int in_offset = i * width * 3;
		for(j = 0; j < width; ++j){
			RGB_To_RGBA((data + out_offset + j*4), (newData + in_offset + j*3));
		}
	}
}

void TargaImage::RGB_To_RGBA(unsigned char *rgba, unsigned char *rgb)
{
	float alpha = rgba[3];
	float alphaScale = alpha / 255.0;
	float val;
	int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = rgb[i] * alphaScale;
	        if (val < 0)
				rgba[i] = 0;
	        else if (val > 255)
				rgba[i] = 255;
	        else
				rgba[i] = val;
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
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

