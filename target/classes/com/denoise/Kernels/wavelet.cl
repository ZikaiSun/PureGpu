
__kernel void haar_wavelet_x (IMAGE_src_TYPE src,
                          IMAGE_dst_TYPE dst,
                          const int nx,
                          const int ny,
                          const int nz,
                          const int inx,
                          const int iny,
                          const int inz
                     )
{
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int z = get_global_id(2);

    const int4 dst_id = (int4)(x, y, z, 0);
    float value;
    // float sqrt_2 = sqrt(2.0F);
    float sqrt_2 = 1;
    int4 src_id1, src_id2;
    int half_x = (int)(nx/2);
    if (x<half_x){
        src_id1 = (int4)(x*2, y, z, 0);
        src_id2 = (int4)(x*2+1, y, z, 0);
        const float value1 = READ_src_IMAGE(src, intsampler, src_id1).x;
        const float value2 = READ_src_IMAGE(src, intsampler, src_id2).x;
        value = (value1+value2)/sqrt_2;
        WRITE_dst_IMAGE(dst, dst_id, CONVERT_dst_PIXEL_TYPE(value));
    }else{
        src_id1 = (int4)((x-half_x)*2, y, z, 0);
        src_id2 = (int4)((x-half_x)*2+1, y, z, 0);
        const float value1 = READ_src_IMAGE(src, intsampler, src_id1).x;
        const float value2 = READ_src_IMAGE(src, intsampler, src_id2).x;
        value = (value1-value2)/sqrt_2;
        WRITE_dst_IMAGE(dst, dst_id, CONVERT_dst_PIXEL_TYPE(value));
    }


}


__kernel void haar_wavelet_y (    IMAGE_src_TYPE src,
                           IMAGE_dst_TYPE dst,
                          const int nx,
                          const int ny,
                          const int nz,
                          const int inx,
                          const int iny,
                          const int inz
                     )
{
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int z = get_global_id(2);

    const int4 dst_id = (int4)(x, y, z, 0);
    float value;
    // float sqrt_2 = sqrt(2.0F);
    float sqrt_2 = 1;
    int4 src_id1, src_id2;

    int half_y = (int)(ny/2);
    if (y<half_y){
        src_id1 = (int4)(x, y*2, z, 0);
        src_id2 = (int4)(x, y*2+1, z, 0);
        const float value1 = READ_src_IMAGE(src, intsampler, src_id1).x;
        const float value2 = READ_src_IMAGE(src, intsampler, src_id2).x;
        value = (value1+value2)/sqrt_2;
        WRITE_dst_IMAGE(dst, dst_id, CONVERT_dst_PIXEL_TYPE(value));
    }else{
        src_id1 = (int4)(x, (y-half_y)*2, z, 0);
        src_id2 = (int4)(x, (y-half_y)*2+1, z, 0);
        const float value1 = READ_src_IMAGE(src, intsampler, src_id1).x;
        const float value2 = READ_src_IMAGE(src, intsampler, src_id2).x;
        value = (value1-value2)/sqrt_2;
        WRITE_dst_IMAGE(dst, dst_id, CONVERT_dst_PIXEL_TYPE(value));
    }
}


__kernel void copy_nxyz (IMAGE_src_TYPE src,
                         IMAGE_dst_TYPE dst,
                          const int src_start_x,
                          const int src_start_y,
                          const int src_start_z,
                          const int dst_start_x,
                          const int dst_start_y,
                          const int dst_start_z
                      )
{
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
    int x = get_global_id(0);
    int y = get_global_id(1);
    int z = get_global_id(2);

    const int4 src_id = (int4)(x+src_start_x,y+src_start_y,z+src_start_z,0);
    const int4 dst_id = (int4)(x+dst_start_x,y+dst_start_y,z+dst_start_z,0);

    const float value = READ_src_IMAGE(src, intsampler, src_id).x;
    WRITE_dst_IMAGE(dst, dst_id, CONVERT_dst_PIXEL_TYPE(value));

}





__kernel void ihaar_wavelet_y (IMAGE_src_TYPE src,
                              IMAGE_dst_TYPE dst,
                              const int nx,
                              const int ny,
                              const int nz,
                              const int inx,
                              const int iny,
                              const int inz
                              )
{
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int z = get_global_id(2);


    const int4 dst_id = (int4)(x, y, z, 0);
    float value;
    float sqrt_2 = 2;

    int half_y = ny/2;
    float y_2 = y/2.0F;

    const int4 src_id1 = (int4)(x, floor(y_2), z, 0);
    const int4 src_id2 = (int4)(x, floor(y_2)+half_y, z, 0);

    const float value1 = READ_src_IMAGE(src, intsampler, src_id1).x;
    const float value2 = READ_src_IMAGE(src, intsampler, src_id2).x;

    if (y%2==0){
        value = (value1+value2)/sqrt_2;
        WRITE_dst_IMAGE(dst, dst_id, CONVERT_dst_PIXEL_TYPE(value));
    }else{
        value = (value1-value2)/sqrt_2;
        WRITE_dst_IMAGE(dst, dst_id, CONVERT_dst_PIXEL_TYPE(value));
    }

}




__kernel void ihaar_wavelet_x (IMAGE_src_TYPE src,
                              IMAGE_dst_TYPE dst,
                              const int nx,
                              const int ny,
                              const int nz,
                              const int inx,
                              const int iny,
                              const int inz
                              )
{
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_NONE | CLK_FILTER_NEAREST;
    const int x = get_global_id(0);
    const int y = get_global_id(1);
    const int z = get_global_id(2);


    const int4 dst_id = (int4)(x, y, z, 0);
    float value;
    float sqrt_2 = 2;

    int half_x = nx/2;
    float x_2 = x/2.0F;

    const int4 src_id1 = (int4)(floor(x_2), y, z, 0);
    const int4 src_id2 = (int4)(floor(x_2)+half_x, y, z, 0);

    const float value1 = READ_src_IMAGE(src, intsampler, src_id1).x;
    const float value2 = READ_src_IMAGE(src, intsampler, src_id2).x;


    if (x%2==0){
        value = (value1+value2)/sqrt_2;
        WRITE_dst_IMAGE(dst, dst_id, CONVERT_dst_PIXEL_TYPE(value));
    }else{
        value = (value1-value2)/sqrt_2;
        WRITE_dst_IMAGE(dst, dst_id, CONVERT_dst_PIXEL_TYPE(value));
    }

}

