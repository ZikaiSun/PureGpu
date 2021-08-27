
__kernel void foldKernel(IMAGE_MN_TYPE MN, int K, int offset) {
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;

    int k1 = get_global_id(0);
    int k2 = get_global_id(1);
    int gid = get_global_id(2);
    float value1 = READ_MN_IMAGE(MN, intsampler, (int4)(k1*K+k2,gid,0,0)).x;
    float value2 = READ_MN_IMAGE(MN, intsampler, (int4)(k1*K+k2,gid+offset,0,0)).x;
    float value = value1 + value2;
    WRITE_MN_IMAGE(MN, (int4)(k1*K+k2,gid,0,0), CONVERT_MN_PIXEL_TYPE(value));

}