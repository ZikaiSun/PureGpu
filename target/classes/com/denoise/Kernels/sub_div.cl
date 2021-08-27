

__kernel void sub_div (IMAGE_src_TYPE src,
                        __global float * alpha,
                        __global float * delta
                      )
{
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
    int x = get_global_id(0);
    int y = get_global_id(1);
    int z = get_global_id(2);

    const int4 src_id = (int4)(x,y,z,0);
    float value = READ_src_IMAGE(src, intsampler, src_id).x;
    value = (value-delta[z])/alpha[z];

    WRITE_src_IMAGE(src, src_id, CONVERT_src_PIXEL_TYPE(value));
    // if(x==1 && y==2 && z==0){
    //     printf("alpha delta %f %f \n", alpha[0], delta[0]);
    // }
}


