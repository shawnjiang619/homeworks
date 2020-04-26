#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853


/***************************** Box filter *******************************
 We want to create a box filter. We will only use square box filters.
 One way to do this is make an image,fill it in with all 1s, and then
 normalize it.That's what we'll do because the normalization function may
 be useful in the future!
 ************************************************************************/
void l1_normalize(image im)
{
    // TODO
    /***********************************************************************
     This function divides each value in an image "im" by the sum of all the
     values in the image and modifies the image in place.
     ************************************************************************/
    // calculate sum
    float sum = 0;
    for (int i = 0; i < im.c; i++) {
        for (int j = 0; j < im.h; j++) {
            for (int k = 0; k < im.w; k++) {
                sum += get_pixel(im, k, j, i);
            }
        }
    }
    if (sum == 0.0) {
        return;
    }
    // divide each value by sum
    for (int i = 0; i < im.c; i++) {
        for (int j = 0; j < im.h; j++) {
            for (int k = 0; k < im.w; k++) {
                float val = get_pixel(im, k, j, i);
                set_pixel(im, k, j, i, val / sum);
            }
        }
    }
}

image make_box_filter(int w)
{
    // TODO
    /***********************************************************************
     This function makes a square filter of size "w x w". Hint:
     Change the "make_image" arguments so that it is a square image of
     width = height = w and number of channels = 1, with all entries equal
     to 1. Then use "l1_normalize" to normalize your filter.
     ************************************************************************/
    // make an image of w * w
    image res = make_image(w, w, 1);
    // fill all entries to 1
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < w; j++) {
            set_pixel(res, i, j, 1, 1);
        }
    }
    // normalize the filter
    l1_normalize(res);
    return res;
}


// a helper function that calculate the appropriate pixel value given a filter
float convHelper (image filter, int channel, image im, int x, int y, int c){
    float val = 0;
    // find position of the upper corner of filter
    int XPos = x - (filter.w-1)/2;
    int YPos = y - (filter.h-1)/2;
    //caluculate the appropriate value
    for (int i = 0; i < filter.w;i++){
        for (int j = 0; j < filter.h;j++){
            // add the product of specific pixel to val
            val += get_pixel(filter, i, j, channel)*get_pixel(im, XPos+i, YPos+j, c);
        }
    }
    return val;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    /***********************************************************************
    im: an image with shape "h x w x c"
    filter: a convolution filter with shape "k1 x k2 x k3". 
    Preserve: an integer, which is usually either 0 or 1.

    - If `filter` and `im` have the same number of channels then it's just a normal 
    convolution. We sum over spatial and channel dimensions and produce a 1 channel image. 
    UNLESS:
        If `preserve` is set to 1 we should produce an image with the same number of 
        channels as the input. This is useful if, for example, we want to run a box 
        filter over an RGB image and get out an RGB image. This means each channel in 
        the image will be filtered by the corresponding channel in the filter.
    If the `filter` only has one channel but `im` has multiple channels we want to
    apply the filter to each of those channels. Then we either sum between channels 
    or not depending on if `preserve` is set.

    Also, `filter` better have either the same number of channels as `im` or have one channel. 
    I check this with an `assert`.
    ************************************************************************/
    assert(filter.c == im.c || filter.c ==1);
    image res = make_image(im.w,im.h,im.c);
    float val;
    // convolve image with given filter
    for(int i = 0; i < im.c; i++){
        for (int j = 0; j < im.w;j++){
            for (int k = 0; k < im.h; k++){
                // check if filter has only one channel
                if (filter.c == 1) {
                    val = convHelper(filter,0, im, j,k,i);
                } else {
                    val = convHelper(filter,i, im, j,k,i);
                }
                // set pixel to calculated value
                set_pixel(res,j,k, i,val);
            }
        }
    }
    // check if we need to preserve image
    if (preserve == 1) {
        return res;
    } else {
        // sum over spatial and channel dimensions and produce a 1 channel image
        image res2 = make_image(im.w,im.h,1);
        for(int i = 0; i < im.w; i++){
            for(int j = 0; j < im.h;j++){
                val = 0;
                for (int k = 0; k < im.c; k++){
                    val += get_pixel(res, i, j ,k);
                }
                set_pixel(res2, i, j ,0, val);
            }
        }
        return res2;
    }
}

image make_highpass_filter()
{
    // TODO
    /***********************************************************************
    Create a 3x3 highpass filter and return it
    ************************************************************************/
    image res = make_image(3,3,1);
    res.data[1] = -1;
    res.data[3] = -1;
    res.data[5] = -1;
    res.data[7] = -1;
    res.data[4] = 4;
    return res;
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
    Create a 3x3 sharpen filter and return it
    ************************************************************************/
    image res = make_image(3,3,1);
    res.data[1] = -1;
    res.data[3] = -1;
    res.data[5] = -1;
    res.data[7] = -1;
    res.data[4] = 5;
    return res;
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
    Create a 3x3 emboss filter and return it
    ************************************************************************/
    image res = make_image(3,3,1);
    res.data[0] = -2;
    res.data[1] = -1;
    res.data[3] = -1;
    res.data[4] = 1;
    res.data[5] = 1;
    res.data[7] = 1;
    res.data[8] = 2;
    return res;
}

// Question 2.2.1: Which of these filters should we use `preserve = 1` when we run our convolution and which ones should we not? Why?
// Answer: TODO
/* we need to use `preserve = 1` for sharpen and emboss filter because we need to preserve colors after applying them.
   On the other hand, we to not use `preserve = 1` for high pass because color is not necessary. It is applied to graytone image.
 */

// Question 2.2.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO
/* Post-processing are necessary for all filters because value of a pixel may become out of range after applying these filters.
 */

image make_gaussian_filter(float sigma)
{
    // TODO
    /***********************************************************************
    sigma: a float number for the Gaussian.

    Create a Gaussian filter with the given sigma. Note that the kernel size 
    is the next highest odd integer from 6x sigma.

    Return the Gaussian filter.
    ************************************************************************/
    // calculate the kernel size
    int size = ceil(sigma * 6);
    if (size % 2 == 0){
        size = size+1;
    }
    image res = make_image(size,size,1);
    // fill the kernel with the probability density function
    for (int i = 0; i < size; i++){
        for (int j = 0; j < size; j++){
            float num = 1.0 / (2 * M_PI * sigma * sigma) * exp(-(i-(size/2) * i-(size/2) + j-(size/2) * j-(size/2)) / (2 * sigma * sigma));
            set_pixel(res,i, j, 0, num);
        }
    }
    //l1_normalize(res);
    return res;
}

image add_image(image a, image b)
{
    // TODO
    /***********************************************************************
    The input images a and image b have the same height, width, and channels.
    Sum the given two images and return the result. 
    The result image should also have the same height, width, and channels as the inputs.
    ************************************************************************/
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image res = make_image(a.w, a.h, a.c);
    for (int i = 0; i < a.w*a.h*a.c; i++){
        res.data[i] = a.data[i] + b.data[i];
    }
    return res;
}

image sub_image(image a, image b)
{
    // TODO
    /***********************************************************************
    The input image a and image b have the same height, width, and channels.
    Subtract the given two images (a - b) and return the result.
    The result image should also have the same height, width, and channels as the inputs.
    ************************************************************************/
    assert(a.w == b.w && a.h == b.h && a.c == b.c);
    image res = make_image(a.w, a.h, a.c);
    for (int i = 0; i < a.w*a.h*a.c; i++){
        res.data[i] = a.data[i] - b.data[i];
    }
    return res;
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
    Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    image res = make_image(3,3,1);
    res.data[0] = -1;
    res.data[2] = 1;
    res.data[3] = -2;
    res.data[5] = 2;
    res.data[6] = -1;
    res.data[8] = 1;
    return res;
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
    Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    image res = make_image(3,3,1);
    res.data[0] = -1;
    res.data[1] = -2;
    res.data[2] = -1;
    res.data[6] = 1;
    res.data[7] = 2;
    res.data[8] = 1;
    return res;
}

image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
    im: the input image is either a "h x w x 3" RGB image or "h x w x 1" grayscale 
    image.
    
    Apply Sobel filter to the given image, get the magnitude and gradient, 
    and return the result. 

    Hint: the "calloc" function can allocate the memory for your output. You can
    assess the first image (magnitute) by calling rst[0] and the second image 
    by calling rst[1]
    ************************************************************************/
    image *rst = calloc(2, sizeof(image));
    
    // make gx and gy filters
    image gx = convolve_image(im, make_gx_filter(),0);
    image gy = convolve_image(im, make_gy_filter(), 0);
    // make magnitute and gradient image
    image magnitute = make_image(im.w, im.h, 1);
    image gradient = make_image(im.w, im.h, 1);
    //compute result
    for (int i = 0; i < gx.w*gx.h;i++){
        magnitute.data[i] = sqrt(gx.data[i] * gx.data[i] + gy.data[i]*gy.data[i]);
        gradient.data[i] = atan2(gy.data[i], gx.data[i]);
    }
    rst[0] = magnitute;
    rst[1] = gradient;
    return rst;
}


void normalize_image(image im)
{
    /***********************************************************************
    Calculate minimum and maximum pixel values. Normalize the image by
    subtracting the minimum and dividing by the max-min difference.
    This is a helper function to visualize the sobel magnitude image better.
    No TODO here :)
    ***********************************************************************/
    int i;
    float min = im.data[0];
    float max = im.data[0];
    for(i = 0; i < im.w*im.h*im.c; ++i){
        if(im.data[i] > max) max = im.data[i];
        if(im.data[i] < min) min = im.data[i];
    }
    for(i = 0; i < im.w*im.h*im.c; ++i){
        im.data[i] = (im.data[i] - min)/(max-min);
    }
}


// EXTRA CREDITS BELOW
int compare_float(const void * a, const void * b)
{
    // This function is provided for your convenience
    float fa = *(const float*) a;
    float fb = *(const float*) b;
    return (fa > fb) - (fa < fb);
}

image apply_median_filter(image im, int kernel_size)
{
    image out = make_image(im.w, im.h, im.c);

    // TODO (EXTRA CREDIT)
    /***********************************************************************
    im is the input image.
    kernel_size is a positive odd number.

    We assume a median filter is a square, with the same height and width.
    The kernel size is always a positive odd number. We use "clamp" padding
    for borders and corners. The output image should have the same width, 
    height, and channels as the input image. You should apply median filter
    to each channel of the input image `im`.

    Hint: use the qsort() function to sort an array. Make use of compare_float() as needed.
    ************************************************************************/

    return out;
}
