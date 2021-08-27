

__kernel void mul_add (IMAGE_src_TYPE src,
                          const float mul,
                          const float add
                      )
{
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
    int x = get_global_id(0);
    int y = get_global_id(1);
    int z = get_global_id(2);

    const int4 id = (int4)(x,y,z,0);
    const float value1 = READ_src_IMAGE(src, intsampler, id).x;

    float value = value1*mul+add;
    
    WRITE_src_IMAGE(src, id, CONVERT_src_PIXEL_TYPE(value));

}


