#include "../image.h"
#include "../utils.h"

#include <string>

using namespace std;

void test_get_pixel()
  {
  Image im = load_image("data/dots.png");
  // Test within image
  TEST(within_eps(0, im.clamped_pixel(0,0,0)));
  TEST(within_eps(1, im.clamped_pixel(1,0,1)));
  TEST(within_eps(0, im.clamped_pixel(2,0,1)));

  // Test padding
  TEST(within_eps(1, im.clamped_pixel(0,3,1)));
  TEST(within_eps(1, im.clamped_pixel(7,8,0)));
  TEST(within_eps(0, im.clamped_pixel(7,8,1)));
  TEST(within_eps(1, im.clamped_pixel(7,8,2)));
  }

void test_set_pixel()
  {
  Image im = load_image("data/dots.png");
  Image d(4,2,3);
  
  d.set_pixel(0,0,0,0); d.set_pixel(0,0,1,0); d.set_pixel(0,0,2,0);
  d.set_pixel(1,0,0,1); d.set_pixel(1,0,1,1); d.set_pixel(1,0,2,1);
  d.set_pixel(2,0,0,1); d.set_pixel(2,0,1,0); d.set_pixel(2,0,2,0);
  d.set_pixel(3,0,0,1); d.set_pixel(3,0,1,1); d.set_pixel(3,0,2,0);
  
  d.set_pixel(0,1,0,0); d.set_pixel(0,1,1,1); d.set_pixel(0,1,2,0);
  d.set_pixel(1,1,0,0); d.set_pixel(1,1,1,1); d.set_pixel(1,1,2,1);
  d.set_pixel(2,1,0,0); d.set_pixel(2,1,1,0); d.set_pixel(2,1,2,1);
  d.set_pixel(3,1,0,1); d.set_pixel(3,1,1,0); d.set_pixel(3,1,2,1);
  
  // Test images are same
  TEST(same_image(im, d));
  }

void test_grayscale()
  {
  Image im = load_image("data/colorbar.png");
  Image gray = rgb_to_grayscale(im);
  Image g = load_image("data/gray.png");
  TEST(same_image(gray, g));
  }

void test_copy()
  {
  Image im = load_image("data/dog.jpg");
  Image c = im;
  TEST(same_image(im, c));
  }

void test_shift()
  {
  Image im = load_image("data/dog.jpg");
  Image c = im;
  shift_image(c, 1, .1);
  TEST(within_eps(im.data[0], c.data[0]));
  TEST(within_eps(im.data[im.w*im.h+13] + .1,  c.data[im.w*im.h + 13]));
  TEST(within_eps(im.data[2*im.w*im.h+72],  c.data[2*im.w*im.h + 72]));
  TEST(within_eps(im.data[im.w*im.h+47] + .1,  c.data[im.w*im.h + 47]));
  shift_image(im, 0, .4);
shift_image(im, 1, .4);
shift_image(im, 2, .4);
im.save_image("data/shift_result");
clamp_image(im);
im.save_image("data/clamp_result");
  }

void test_scale()
{
    Image im = load_image("data/dog.jpg");
    Image s = im;
    rgb_to_hsv(im);
    rgb_to_hsv(s);
    scale_image(s, 1, 2);
    int ch_size = im.w*im.h;
    TEST(within_eps(im.data[ch_size]*2, s.data[ch_size]));
    TEST(within_eps(im.data[ch_size+72]*2,  s.data[ch_size + 72]));
    TEST(within_eps(im.data[2*ch_size+72],  s.data[2*ch_size + 72]));
    TEST(within_eps(im.data[47],  s.data[47]));
    
 im = load_image("data/dog.jpg");
rgb_to_hsv(im);
scale_image(im, 1, 2);
clamp_image(im);
hsv_to_rgb(im);
im.save_image("data/dog_scale_saturated");

}

void test_rgb_to_hsv()
  {
  Image im = load_image("data/dog.jpg");
  rgb_to_hsv(im);
  Image hsv = load_image("data/dog.hsv.png");
  TEST(same_image(im, hsv));
  }

void test_hsv_to_rgb()
  {
  Image im = load_image("data/dog.jpg");
  Image c = im;
  rgb_to_hsv(im);
  hsv_to_rgb(im);
  TEST(same_image(im, c));
  Image im2 = load_image("data/dog.jpg");
rgb_to_hsv(im2);
shift_image(im2, 1, .2);
clamp_image(im2);
hsv_to_rgb(im2);
im2.save_image("data/colorspace_result");
  }

void test_rgb2lch2rgb()
  {
  Image im = load_image("data/dog.jpg");
  Image c = im;
  
  rgb_to_lch(im);
  lch_to_rgb(im);
  TEST(same_image(im, c));
  }


void run_tests()
  {
  test_get_pixel();
  test_set_pixel();
  test_copy();
  test_shift();
  test_scale();
  test_grayscale();
  test_rgb_to_hsv();
  test_hsv_to_rgb();
//  test_rgb2lch2rgb();
  printf("%d tests, %d passed, %d failed\n", tests_total, tests_total-tests_fail, tests_fail);
  }

int main(int argc, char **argv)
  {
  // Image manipulation for fun testing.
  
  Image im2 = load_image("data/dog.jpg");
  for (int i=0; i<im2.w; i++)
      for (int j=0; j<im2.h; j++)
          im2(i, j, 0) = 0;
  im2.save_image("data/pixel_modifying_output");
  
  // Running example tests
  
  run_tests();
  
  return 0;
  }
