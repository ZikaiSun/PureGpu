

__kernel void get_parent(IMAGE_low_TYPE low, IMAGE_dst_TYPE dst,
                        int offset_flag,
                      int nx, int ny, int nz)
{
  const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;

    int x = get_global_id(0);
    int y = get_global_id(1);
    int z = get_global_id(2);


  float value;
  if (offset_flag == 1 ){
      // int ym1 = (y-1)>=0?(y-1):y;
      // int yp1 = (y+1)<ny?(y+1):y;
      const int4 src_id1 = (int4)(x,y-1,z,0);
      const int4 src_id2 = (int4)(x,y+1,z,0);
      value = READ_low_IMAGE(low, intsampler, src_id1).x-READ_low_IMAGE(low, intsampler, src_id2).x;
  }
  else if (offset_flag == 2){
      // int xm1 = (x-1)>=0?(x-1):x;
      // int xp1 = (x+1)<nx?(x+1):x;
      const int4 src_id1 = (int4)(x-1,y,z,0);
      const int4 src_id2 = (int4)(x+1,y,z,0);
      value = READ_low_IMAGE(low, intsampler, src_id1).x-READ_low_IMAGE(low, intsampler, src_id2).x;
  }
  else if (offset_flag == 3){
      // int xm1 = (x-1)>=0?(x-1):x;
      // int xp1 = (x+1)<nx?(x+1):x;
      // int ym1 = (y-1)>=0?(y-1):y;
      // int yp1 = (y+1)<ny?(y+1):y;
      const int4 src_id1 = (int4)(x-1,y-1,z,0);
      const int4 src_id2 = (int4)(x+1,y-1,z,0);
      const int4 src_id3 = (int4)(x-1,y+1,z,0);
      const int4 src_id4 = (int4)(x+1,y+1,z,0);
      value = READ_low_IMAGE(low, intsampler, src_id1).x-READ_low_IMAGE(low, intsampler, src_id2).x-READ_low_IMAGE(low, intsampler, src_id3).x+READ_low_IMAGE(low, intsampler, src_id4).x;
  }
  else{
    value = 0;
  }

  const int4 dst_id = (int4)(x,y,z,0);
  WRITE_dst_IMAGE(dst, dst_id, CONVERT_dst_PIXEL_TYPE(value));

}


// __kernel void get_parent(IMAGE_low_TYPE low, IMAGE_dst_TYPE dst,
//                         int offset_flag,
//                       int nx, int ny, int nz)
// {
//   const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
//
//     int x = get_global_id(0);
//     int y = get_global_id(1);
//     int z = get_global_id(2);
//
//
//   float value;
//   if (offset_flag == 1 ){
//       int ym1 = (y-1)>=0?(y-1):y;
//       int yp1 = (y+1)<ny?(y+1):y;
//       const int4 src_id1 = (int4)(x,ym1,z,0);
//       const int4 src_id2 = (int4)(x,yp1,z,0);
//       value = READ_low_IMAGE(low, intsampler, src_id1).x-READ_low_IMAGE(low, intsampler, src_id2).x;
//   }
//   else if (offset_flag == 2){
//       int xm1 = (x-1)>=0?(x-1):x;
//       int xp1 = (x+1)<nx?(x+1):x;
//       const int4 src_id1 = (int4)(xm1,y,z,0);
//       const int4 src_id2 = (int4)(xp1,y,z,0);
//       value = READ_low_IMAGE(low, intsampler, src_id1).x-READ_low_IMAGE(low, intsampler, src_id2).x;
//   }
//   else if (offset_flag == 3){
//       int xm1 = (x-1)>=0?(x-1):x;
//       int xp1 = (x+1)<nx?(x+1):x;
//       int ym1 = (y-1)>=0?(y-1):y;
//       int yp1 = (y+1)<ny?(y+1):y;
//       const int4 src_id1 = (int4)(xm1,ym1,z,0);
//       const int4 src_id2 = (int4)(xp1,ym1,z,0);
//       const int4 src_id3 = (int4)(xm1,yp1,z,0);
//       const int4 src_id4 = (int4)(xp1,yp1,z,0);
//       value = READ_low_IMAGE(low, intsampler, src_id1).x-READ_low_IMAGE(low, intsampler, src_id2).x-READ_low_IMAGE(low, intsampler, src_id3).x+READ_low_IMAGE(low, intsampler, src_id4).x;
//   }
//   else{
//     value = 0;
//   }
//
//   const int4 dst_id = (int4)(x,y,z,0);
//   WRITE_dst_IMAGE(dst, dst_id, CONVERT_dst_PIXEL_TYPE(value));
//
// }
