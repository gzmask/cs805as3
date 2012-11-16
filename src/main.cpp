#include <iostream>
#include "util.h"

/* definition of window on the image plane in the camera coordinates */
/* They are used in mapping (j, i) in the screen coordinates into */
/* (x, y) on the image plane in the camera coordinates */
/* The window size used here simulates the 35 mm film. */
double xmin = 0.0175;
double ymin = -0.0175;
double xmax = -0.0175;
double ymax = 0.0175;

/* camera position */
Point VRP = {128.0, 64.0, 250.0};
Vector VPN = {-64.0, 0.0, -186.0};
Vector VUP = {0.0, 1.0, 0.0};

double focal = 0.05;	/* focal length simulating 50 mm lens */

Vector Light = {0.577, -0.577, -0.577}; /* light direction */
double Ip = 255.0; /* intensity of the point light source */

/* Transformation from the world to the camera coordinates */
Matrix Mwc = get_M(VRP, VPN, VUP);
Matrix Rwc = get_R(VRP, VPN, VUP);
Matrix Twc = get_T(VRP);
/* Transformation from the camera to the world coordinates */
Matrix Mcw = get_Mi(VRP, VPN, VUP);
Matrix Rcw = get_Ri(VRP, VPN, VUP);
Matrix Tcw = get_Ti(VRP);

int main () {

  //main program for volume rendering
  Volume* ct = new Volume;
  read_from_file("smallHead.den", ct);
  //print_ct_volume(ct);
  Volume* color = new Volume;
  compute_shading_volume(ct, color);
  //print_ct_volume(color);
  ImagePanel* img = new ImagePanel;
  init_img_panel(img);
  foreach_pixel_exec(img, volume_ray_tracing, ct, color);
  save_to_file(img, "output.raw"); //save result to binary file
  delete img;
  delete ct;
  delete color;

  return 0;
}
