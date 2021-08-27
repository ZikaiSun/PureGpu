
__kernel void getC(IMAGE_tmp_TYPE tmp, int K, int ny, int offset) {
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;

    int k0 = get_global_id(0);
    int gid = get_global_id(1);
    int x1 = gid%ny;
    int y1 = gid/ny;
    int x2 = (gid+offset)%ny;
    int y2 = (gid+offset)/ny;

    float value1 = READ_tmp_IMAGE(tmp, intsampler, (int4)(k0,x1,y1,0)).x;
    float value2 = READ_tmp_IMAGE(tmp, intsampler, (int4)(k0,x2,y2,0)).x;
    float value = value1 + value2;
    WRITE_tmp_IMAGE(tmp, (int4)(k0,x1,y1,0), CONVERT_tmp_PIXEL_TYPE(value));

}