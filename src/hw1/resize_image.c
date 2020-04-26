#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"


/***********************************************************************
  We've been talking a lot about resizing and interpolation in class,
  now's your time to do it!
  In order to make grading easier, please only edit the files we mention to submit.
  You will submit the resize_image.c file on Canvas.
************************************************************************/


/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

float nn_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs nearest-neighbor interpolation on image "im"
      given a floating column value "x", row value "y" and integer channel "c",
      it interpolates and returns the interpolated value.
      Remember to use the closest "int", not just type-cast because in C that
      will truncate towards zero.
    ************************************************************************/
    // find nearest neighbor
    int lx = round(x);
    int ly = round(y);
    // get its pixel
    float res = get_pixel(im, lx, ly, c);
    return res;
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix the return line)
    /***********************************************************************
      This function uses nearest-neighbor interpolation on image "im" to
      construct a new image of size "w x h". Hint:
      - Create a new image that is "w x h" and the same number of channels as "im"
      - Loop over the pixels and map back to the old coordinates.
      - Use nearest-neighbor interpolate to fill in the image.
    ************************************************************************/
    // create image with specific w and h
    image res = make_image(w, h, im.c);
    // caluculate the scale of new and old image
    float w_a = (float)im.w / w;
    float w_b = -0.5 + 0.5 * w_a;
    float h_a = (float)im.h / h;
    float h_b = -0.5 + 0.5 * h_a;
    for (int k = 0; k < im.c; k++) {
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                float x = w_a * i + w_b;
                float y = h_a * j + h_b;
                // calculate new pixel value
                float val = nn_interpolate(im, x, y, k);
                set_pixel(res, i, j, k, val);
            }
        }
    }
    return res;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs bilinear interpolation on image "im" given
      a floating column value "x", row value "y" and integer channel "c".
      It interpolates and returns the interpolated value.
    ************************************************************************/
    // find distances
    float dx = x - floor(x);
    float dy = y - floor(y);
    // calculate four nearby points
    float up_left = get_pixel(im, floor(x), floor(y), c);
    float up_right = get_pixel(im, floor(x)+1, floor(y), c);
    float down_left = get_pixel(im, floor(x), floor(y)+1, c);
    float down_right = get_pixel(im, floor(x)+1, floor(y)+1, c);
    float res =   up_left*(1-dx)*(1-dy) + up_right*dx*(1-dy) + down_left*(1-dx)*dy + down_right*dx*dy;
    return res;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    /***********************************************************************
      This function uses bilinear interpolation on image "im" to construct
      a new image of size "w x h". Hint:
      - Create a new image that is "w x h" and the same number of channels as "im".
      - Loop over the pixels and map back to the old coordinates.
      - Use bilinear interpolate to fill in the image.
    ************************************************************************/
    // caluculate the scale of new and old image
    float w_a = (float)im.w / w;
    float w_b = -0.5 + 0.5 * w_a;
    float h_a = (float)im.h / h;
    float h_b = -0.5 + 0.5 * h_a;
    // create new image
    image res = make_image(w, h, im.c);
    for (int k = 0; k < im.c; k++) {
        for (int j = 0; j < h; j++) {
            for (int i = 0; i < w; i++) {
                float x = w_a * i + w_b;
                float y = h_a * j + h_b;
                // calculate new pixel value
                float val = bilinear_interpolate(im, x, y, k);
                set_pixel(res, i, j, k, val);
            }
        }
    }
    return res;
}
