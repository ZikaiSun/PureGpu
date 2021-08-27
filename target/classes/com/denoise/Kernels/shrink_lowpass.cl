#pragma OPENCL EXTENSION cl_khr_fp64 : enable


__kernel void shrink(IMAGE_low_TYPE low, IMAGE_F_TYPE F, IMAGE_tmp_TYPE tmp, IMAGE_sigma2_TYPE sigma2,
                     int currentFrame, int nx, int ny, int nz)
{
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;

    int x = get_global_id(0);
    int y = get_global_id(1);

    int N = nx * ny;
    float Yl, yl, ylm;
    float EPS = 1e-9;
    float Ec[20] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
    Ec[currentFrame]=1;


    int4 src_id = (int4)(x,y,currentFrame,0);
    Yl = READ_low_IMAGE(low, intsampler, src_id).x;
    float sigma2_currentframe = READ_sigma2_IMAGE(sigma2, intsampler, (int4)(currentFrame,0,0,0)).x;


    int4 ijc_id;
    int4 cxy_id;
    float value;

    for (int c = 0; c < nz; c++) {
        ijc_id = (int4)(x,y,c,0);
        yl = READ_low_IMAGE(low, intsampler, ijc_id).x;
        ylm = yl - Ec[c];

        cxy_id = (int4)(c,x,y,0);

        WRITE_F_IMAGE(F, cxy_id, CONVERT_F_PIXEL_TYPE(yl));
        value = ((Yl * ylm - sigma2_currentframe * Ec[c]) / N);
        WRITE_tmp_IMAGE(tmp, cxy_id, CONVERT_tmp_PIXEL_TYPE(value));
    }

}

