

__kernel void sum_up (IMAGE_src_TYPE src,
                         IMAGE_dst_TYPE dst
                      )
{
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
    int x = get_global_id(0);
    int y = get_global_id(1);
    int z = get_global_id(2);

    const int4 id = (int4)(x,y,z,0);
    const float value1 = READ_src_IMAGE(src, intsampler, id).x;
    const float value2 = READ_dst_IMAGE(dst, intsampler, id).x;
    float value = value1 + value2;

    WRITE_dst_IMAGE(dst, id, CONVERT_dst_PIXEL_TYPE(value));

}


