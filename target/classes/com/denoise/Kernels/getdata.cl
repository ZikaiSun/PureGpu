

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


