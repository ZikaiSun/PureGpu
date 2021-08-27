

__kernel void shrink_high(IMAGE_high_TYPE high, IMAGE_low_TYPE low, IMAGE_parent_TYPE parent, IMAGE_F_TYPE F, IMAGE_tmp_TYPE tmp,
                            IMAGE_sigma2_TYPE sigma2,
                        int currentFrame,
                      int nx, int ny, int nz)
{
  const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;

    int x = get_global_id(0);
    int y = get_global_id(1);
    // int x_debug = 31;
    // int y_debug = 5;


    float parent_z[20];
    float smoothed_z[20];

    float sigma=1.0F;
    float Ec[20] = {0}; // nz
    Ec[currentFrame]=1;

    int4 src_id;

    for (int z=0; z<nz; z++){
      src_id = (int4)(x,y,z,0);
      parent_z[z] = READ_parent_IMAGE(parent, intsampler, src_id).x;

        // if (x==x_debug && y==y_debug){
        //     printf("new parent_z[z]: %f \n", parent_z[z]);
        // }

      // smooth
      const int   N = 1;
      const int   c = (N-1)/2;

      float smoothed = 0;
      float sum = 0;
      int x_loc, y_loc;

      for (int a=-c; a<=c; a++){
          for (int b=-c; b<=c; b++){
              float temp = exp(-((float)(a*a+b*b) / (2*sigma*sigma)));
              // x_loc = (x+a)<0?0:( (x+a)>(nx-1)?(nx-1):(x+a) );
              // y_loc = (y+b)<0?0:( (y+b)>(ny-1)?(ny-1):(y+b) );
              x_loc = (x+a);
              y_loc = (y+b);
              src_id = (int4)(x_loc,y_loc,z,0);
              float value_pos = READ_parent_IMAGE(parent, intsampler, src_id).x;

              if (value_pos < 0){
                  value_pos = -1 * value_pos;
              }
              smoothed += temp * value_pos;
              sum += temp;
          }
      }
      smoothed /= sum;
      smoothed_z[z] = smoothed;
    }



    //shrink
    int N = nx * ny, k;
    float yh, yp, yhm, yhp;
    float EPS = 1e-9;
    int4 ijc_id;
    float value;

    src_id = (int4)(x,y,currentFrame,0);
    float Yh = READ_high_IMAGE(high, intsampler, src_id).x;
    float Yl = READ_low_IMAGE(low, intsampler, src_id).x;

    float Ph_0, Ph_1, Ph_2;
    float Dd_0, Dd_1, Dd_2;
    float Pl_0, Pl_1, Pl_2;

    float pixel;
    float yh2 = 0.0;
    float yhm2 = 0.0;
    float yhp2 = 0.0;
    float yp2 = 0.0;
    float var = 0.0;
    float varm = 0.0;
    float T2 = 0.0;
    float T2m = 0.0;
    float ypm2;
    for (int c = 0; c < nz; c++) {
        ijc_id = (int4)(x,y,c,0);
        pixel = READ_high_IMAGE(high, intsampler, ijc_id).x;
        yh2 = yh2 + pixel * pixel;
        yhm2 = yhm2 + (pixel - Ec[c]) * (pixel - Ec[c]);
        yhp2 = yhp2 + (pixel + Ec[c]) * (pixel + Ec[c]);

        pixel = smoothed_z[c];
        yp2 = yp2 + pixel * pixel;

        pixel = READ_low_IMAGE(low, intsampler, ijc_id).x;
        var = var + pixel;
        varm = varm + pixel - Ec[c];

        float sigma2_c = READ_sigma2_IMAGE(sigma2, intsampler, (int4)(c,0,0,0)).x;
        T2 = T2 + sigma2_c;
        T2m = T2m + sigma2_c;
    }
    T2 = (T2 + fabs(var) + EPS) / sqrt((float)nz);
    T2m = (T2m + fabs(varm) + EPS) / sqrt((float)nz);
    yh2 = exp(-(1.0 / 12.0) * yh2 / T2);
    yhm2 = exp(-(1.0 / 12.0) * yhm2 / T2m);
    yhp2 = exp(-(1.0 / 12.0) * yhp2 / T2m);
    ypm2 = exp(-(1.0 / 12.0) * yp2 / T2m);
    yp2 = exp(-(1.0 / 12.0) * yp2 / T2);

    k = 0;
    int k_max = nz*3;
    float Sh_0, Sl_0, DD_0;
    float sigma2_currentframe = READ_sigma2_IMAGE(sigma2, intsampler, (int4)(currentFrame,0,0,0)).x;

    for (int c = 0; c < nz; c++) {
        ijc_id = (int4)(x,y,c,0);
        yh = READ_high_IMAGE(high, intsampler, ijc_id).x;
        yp = parent_z[c];


        yhm = yh - Ec[c];
        yhp = yh + Ec[c];
        // ---------------k=0---------------
        WRITE_F_IMAGE(F, (int4)(k,x,y,0), CONVERT_F_PIXEL_TYPE(yh));

        Ph_0 = yhm + yhp;
        Pl_0 = yhm - yhp;
        Dd_0 = 2.0 * Ec[c];
        Sh_0 = Yh * Ph_0;
        Sl_0 = Yl * Pl_0;
        DD_0 = Dd_0;

        value = (0.5 * (Sh_0 + Sl_0 - sigma2_currentframe * DD_0) / N);
        WRITE_tmp_IMAGE(tmp, (int4)(k,x,y,0), CONVERT_tmp_PIXEL_TYPE(value));

        value = yh * yp2;
        WRITE_F_IMAGE(F, (int4)(k+k_max,x,y,0), CONVERT_F_PIXEL_TYPE(value));

        Sh_0 = Yh * Ph_0 * ypm2;
        Sl_0 = Yl * Pl_0 * ypm2;
        DD_0 = Dd_0 * ypm2 + 2.0 * Pl_0 / (12.0 * T2m * T2m) * ypm2 *  ((varm > 0) ? 1 : ((varm < 0) ? -1 : 0));
        value = (0.5 * (Yh*Ph_0*ypm2 + Yl*Pl_0*ypm2 - sigma2_currentframe * DD_0) / N);
        WRITE_tmp_IMAGE(tmp, (int4)(k+k_max,x,y,0), CONVERT_tmp_PIXEL_TYPE(value));

        // ---------------k=1---------------
        k++;
        value = yh * yh2;
        WRITE_F_IMAGE(F, (int4)(k,x,y,0), CONVERT_F_PIXEL_TYPE(value));

        Ph_1 = yhm * yhm2 + yhp * yhp2;
        Pl_1 = yhm * yhm2 - yhp * yhp2;
        Dd_1 = (Ec[c] - yhm * (Yh - 1.0) / (6.0 * T2m)) * yhm2
                + (Ec[c] - yhp * (Yh + 1.0) / (6.0 * T2m)) * yhp2
                + (yhm * yhm2 - yhp * yhp2) / (12.0 * T2m * T2m) *  ((varm > 0) ? 1 : ((varm < 0) ? -1 : 0));
        Sh_0 = Yh * Ph_1;
        Sl_0 = Yl * Pl_1;
        DD_0 = Dd_1;
        value = (0.5 * (Sh_0 + Sl_0 - sigma2_currentframe * DD_0) / N);
        WRITE_tmp_IMAGE(tmp, (int4)(k,x,y,0), CONVERT_tmp_PIXEL_TYPE(value));

        value = yh * yh2 * yp2;
        WRITE_F_IMAGE(F, (int4)(k+k_max,x,y,0), CONVERT_F_PIXEL_TYPE(value));


        Sh_0 = Yh * Ph_1 * ypm2;
        Sl_0 = Yl * Pl_1 * ypm2;
        DD_0 = Dd_1 * ypm2 + 2.0 * Pl_1 / (12.0 * T2m * T2m) * ypm2 *  ((varm > 0) ? 1 : ((varm < 0) ? -1 : 0));
        value = (0.5 * (Sh_0 + Sl_0 - sigma2_currentframe * DD_0) / N);
        WRITE_tmp_IMAGE(tmp, (int4)(k+k_max,x,y,0), CONVERT_tmp_PIXEL_TYPE(value));

        // ---------------k=2---------------
        k++;
        value = yp;
        WRITE_F_IMAGE(F, (int4)(k,x,y,0), CONVERT_F_PIXEL_TYPE(value));

        Ph_2 = 2.0 * yp;
        Pl_2 = 0;
        Dd_2 = 0;
        Sh_0 = Yh * Ph_2;
        Sl_0 = 0;
        DD_0 = 0;
        value = (0.5 * (Sh_0 + Sl_0 - sigma2_currentframe * DD_0) / N);
        // if (x==x_debug && y==y_debug){
        //     printf("new Yh yp sh sl sigma dd n: %f %f %f %f %f %f %f %f \n", Yh, yp, Sh_0, Sl_0, sigma2_currentframe, DD_0, N, value);
        //     printf("new c_loop %d %d %d %f \n", k, x, y, value);
        // }

        WRITE_tmp_IMAGE(tmp, (int4)(k,x,y,0), CONVERT_tmp_PIXEL_TYPE(value));

        value = yp * yp2;
        WRITE_F_IMAGE(F, (int4)(k+k_max,x,y,0), CONVERT_F_PIXEL_TYPE(value));

        Sh_0 = Yh * Ph_2 * ypm2;
        Sl_0 = Yl * Pl_2 * ypm2;
        DD_0 = Dd_2 * ypm2 + 2.0 * Pl_2 / (12.0 * T2m * T2m) * ypm2 *  ((varm > 0) ? 1 : ((varm < 0) ? -1 : 0));
        value = (0.5 * (Sh_0 + Sl_0 - sigma2_currentframe * DD_0) / N);
        WRITE_tmp_IMAGE(tmp, (int4)(k+k_max,x,y,0), CONVERT_tmp_PIXEL_TYPE(value));

        // ---------------k=3---------------
        k++;
    }

    // if (x==x_debug && y==y_debug){
    //     printf("new end: %d %d %f \n", x, y, READ_tmp_IMAGE(tmp, intsampler, (int4)(2,x,y,0)).x);
    // }

}
