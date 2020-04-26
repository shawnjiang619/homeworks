#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

float get_pixel(image im, int x, int y, int c)
{
    // TODO Fill this in
    // convert input using clamp padding
    if(x >= im.w) x = im.w - 1;
    if(y >= im.h) y = im.h - 1;
    if(x < 0) x = 0;
    if(y < 0) y = 0;
    return im.data[c*im.w*im.h + y*im.w + x];
}

void set_pixel(image im, int x, int y, int c, float v)
{
    // TODO Fill this in
    // check if all parameters are valid
    if (x < 0 || y < 0 || x >= im.w || y >= im.h || c < 0 || c>=im.c) return;
    im.data[x + im.w*y + im.w*im.h*c] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    // TODO Fill this in
    memcpy(copy.data, im.data, im.w*im.h*im.c*sizeof(float));
    return copy;
}

image rgb_to_grayscale(image im)
{
    assert(im.c == 3);
    image gray = make_image(im.w, im.h, 1);
    // TODO Fill this in
    for (int i = 0; i < im.w; i++) {
        for (int j = 0; j < im.h; j++){
            // calculate luma value
            gray.data[i*im.h+j] = 0.299 * im.data[i*im.h+j] + 0.587 * im.data[im.w*im.h +i*im.h+j] + 0.114 * im.data[2*im.w*im.h + i*im.h+j];
        }
    }
    return gray;
}

void shift_image(image im, int c, float v)
{
    // TODO Fill this in
    for(int i = 0; i < im.w*im.h; i++){
        im.data[i + c*im.w*im.h] += v;
    }
}


void scale_image(image im, int c, float v)
{
    for(int i = 0; i < im.w*im.h; i++){
        im.data[i + c*im.w*im.h] *= v;
    }
}

void clamp_image(image im)
{
    // TODO Fill this in
    for(int i = 0; i < im.c*im.w*im.h; i++){
        // make sure the pixel values in the image stay between 0 and 1
        im.data[i] = (im.data[i] < 0) ? 0 : ((im.data[i] > 1) ? 1 : im.data[i]);
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

void rgb_to_hsv(image im)
{
    // TODO Fill this in
    float r, g, b;
    float h, s, v;
    for(int i = 0; i < im.h; i++){
        for(int j = 0; j < im.w; j++){
            // get rgb value
            r = get_pixel(im, j , i, 0);
            g = get_pixel(im, j , i, 1);
            b = get_pixel(im, j , i, 2);
            v = three_way_max(r,g,b);
            //caluculate how far apart the min and max are:
            float m = three_way_min(r,g,b);
            float C = v - m;
            if(v == 0){
                // we do not need extra calculation
                s = 0;
                h = 0;
            }else{
                // caluculate h based on formula
                s = C/v;
                if (C == 0) {
                    h = 0;
                } else if(v == r){
                    h = (g - b) / C;
                } else if (v == g) {
                    h = (b - r) / C + 2;
                } else {
                    h = (r - g) / C  + 4;
                }
                if (h < 0) {
                    h = h/6.0+1;
                } else {
                    h = h/6.0;
                }
            }
            
            // set hsv value respectively
            set_pixel(im, j, i, 0, h);
            set_pixel(im, j, i, 1, s);
            set_pixel(im, j, i, 2, v);
        }
    }
}

void hsv_to_rgb(image im)
{
    // TODO Fill this in
    float r, g, b;
    float h, s, v;
    float f, p, q, t;
    for(int i = 0; i < im.h; i++){
        for(int j = 0; j < im.w; j++){
            // fetch original hsv values
            h = get_pixel(im, j , i, 0);
            s = get_pixel(im, j , i, 1);
            v = get_pixel(im, j , i, 2);
            h = h * 6;
            if (s == 0) {
                r = g = b = v;
            } else {
                int hi = floor(h);
                f = h - hi;
                p = v*(1-s);
                q = v*(1-f*s);
                t = v*(1-(1-f)*s);
                // find rgb values by case
                if(hi == 0){
                    r = v; g = t; b = p;
                } else if(hi == 1){
                    r = q; g = v; b = p;
                } else if(hi == 2){
                    r = p; g = v; b = t;
                } else if(hi == 3){
                    r = p; g = q; b = v;
                } else if(hi == 4){
                    r = t; g = p; b = v;
                } else {
                    r = v; g = p; b = q;
                }
            }
            // set rgb value respectively
            set_pixel(im, j, i, 0, r);
            set_pixel(im, j, i, 1, g);
            set_pixel(im, j, i, 2, b);
        }
    }
}
