
__kernel void getOutput (IMAGE_F_TYPE F, IMAGE_A_TYPE A, IMAGE_out_TYPE out,  IMAGE_indices_TYPE indices, int K)
{
    const sampler_t intsampler  = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_REPEAT | CLK_FILTER_NEAREST;
    const int i = get_global_id(0);
    const int j = get_global_id(1);

    float F_kij, A_k;
    float value = 0.0;
    for (int k = 0; k < K; k++){
        float indices_k = READ_indices_IMAGE(indices, intsampler, (int4)(k,0,0,0)).x;
        int indices_k_int = (int)indices_k;
        if (indices_k_int != -1){
            F_kij = READ_F_IMAGE(F, intsampler, (int4)(indices_k_int,i,j,0)).x;
            A_k = READ_A_IMAGE(A, intsampler, (int4)(k,0,0,0)).x;
            value += F_kij * A_k;
        }
    }
    WRITE_out_IMAGE(out, (int4)(i,j,0,0), CONVERT_out_PIXEL_TYPE(value));

}

