

__kernel void getMN(IMAGE_F_TYPE F, IMAGE_MN_TYPE MN,
                       int K,
                       int nx,
                       int ny
                       )
{
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;


    int k = get_global_id(0);
    int i = get_global_id(1);
    int j = get_global_id(2);

    int k1 = k%K;
    int k2 = k/K;

    float value1 = READ_F_IMAGE(F, intsampler, (int4)(k1,i,j,0)).x;
    float value2 = READ_F_IMAGE(F, intsampler, (int4)(k2,i,j,0)).x;
    float value = value1 * value2;
    WRITE_MN_IMAGE(MN, (int4)(k1*K+k2,i*ny+j,0,0), CONVERT_MN_PIXEL_TYPE(value));

}
