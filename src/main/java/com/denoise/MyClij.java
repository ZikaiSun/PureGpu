package com.denoise;

//====================================================================
//
// Organization:
// Author: Zikai Sun
// Email: zksun@link.cuhk.edu.hk
// Image and Video Processing Lab, The Chinese University of Hong Kong
// Apr 2021
//
// Original Informations:
// Project: PureDenoise
// Package: denoise
// Class  : Denoising
//
// Organization:
// Florian Luisier
// Biomedical Imaging Group (BIG)
// Ecole Polytechnique Fédérale de Lausanne (EPFL)
// Lausanne, Switzerland
//
// Information:
// http://bigwww.epfl.ch/algorithms/denoise/
//
// References:
// [1]  F. Luisier, C. Vonesch, T. Blu, M. Unser, "Fast Interscale Wavelet
//      Denoising of Poisson-corrupted Images", Signal Processing, vol. 90,
//      no. 2, pp. 415-427, February 2010.
// [2]  F. Luisier, "The SURE-LET Approach to Image Denoising," Swiss Federal
//      Institute of Technology Lausanne, EPFL Thesis no. 4566 (2010), 232 p.,
//      January 8, 2010.
// [3]  F. Luisier, C. Vonesch, T. Blu, M. Unser, "Fast Haar-Wavelet Denoising
//      of Multidimensional Fluorescence Microscopy Data", Proceedings of the
//      Sixth IEEE International Symposium on Biomedical Imaging: From Nano to
//      Macro (ISBI'09)}, Boston MA, USA, June 28-July 1, 2009, pp. 310-313.
//
// Conditions of use:
// You'll be free to use this software for research purposes, but you
// should not redistribute it without our consent. In addition, we
// expect you to include a citation or acknowledgement whenever
// you present or publish results that are based on it.
//
//====================================================================

import ij.ImagePlus;
import ij.ImageStack;
import imageware.Builder;
import imageware.ImageWare;
import Jama.LUDecomposition;
import Jama.Matrix;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij.clearcl.interfaces.ClearCLImageInterface;
import net.haesleinhuepf.clij2.CLIJ2;

import java.util.HashMap;

import static net.haesleinhuepf.clij.utilities.CLIJUtilities.assertDifferent;
import ij.IJ;


public class MyClij {

    public static void extend_dim(CLIJ2 clij2, ClearCLImageInterface src, ClearCLImageInterface dst,
                                  int src_start_x, int src_start_y, int src_start_z,
                                  int dst_start_x, int dst_start_y, int dst_start_z,
                                  int nx, int ny, int nz) {
        HashMap<String, Object> parameters1 = new HashMap<>();
        parameters1.put("src", src);
        parameters1.put("dst", dst);
        parameters1.put("src_start_x", src_start_x);
        parameters1.put("src_start_y", src_start_y);
        parameters1.put("src_start_z", src_start_z);
        parameters1.put("dst_start_x", dst_start_x);
        parameters1.put("dst_start_y", dst_start_y);
        parameters1.put("dst_start_z", dst_start_z);
        long globalWorkSize1[] = {nx, ny, nz};
        clij2.execute(MyClij.class, "Kernels/extend_dim.cl", "extend_dim",
                globalWorkSize1, globalWorkSize1, parameters1);
    }

    public static void sub_div(CLIJ2 clij2, ClearCLImageInterface src,
                               int nx, int ny, int nz,
                               double[] alpha, double[] delta, int[] index) {

        int len_stack = index[1] - index[0] + 1;

        float[] alpha_arr = new float[len_stack];
        float[] delta_arr = new float[len_stack];
        for (int i=0; i<len_stack; i++){
            alpha_arr[i] = (float) alpha[index[0]+i];
            delta_arr[i] = (float) delta[index[0]+i];
        }
        ClearCLBuffer alpha_clij = clij2.pushMatXYZ(alpha_arr);
        ClearCLBuffer delta_clij = clij2.pushMatXYZ(delta_arr);

        HashMap<String, Object> parameters1 = new HashMap<>();
        parameters1.put("src", src);
        parameters1.put("alpha", alpha_clij);
        parameters1.put("delta", delta_clij);

        long globalWorkSize1[] = {nx, ny, nz};
        clij2.execute(MyClij.class, "Kernels/sub_div.cl", "sub_div",
                globalWorkSize1, globalWorkSize1, parameters1);
    }

    public static void get_data_periodic(CLIJ2 clij2, ClearCLImageInterface src, ClearCLImageInterface dst,
                                         int src_start_x, int src_start_y, int src_start_z,
                                         int dst_start_x, int dst_start_y, int dst_start_z,
                                         int nx, int ny, int nz) {
        HashMap<String, Object> parameters1 = new HashMap<>();
        parameters1.put("src", src);
        parameters1.put("dst", dst);
        parameters1.put("src_start_x", src_start_x);
        parameters1.put("src_start_y", src_start_y);
        parameters1.put("src_start_z", src_start_z);
        parameters1.put("dst_start_x", dst_start_x);
        parameters1.put("dst_start_y", dst_start_y);
        parameters1.put("dst_start_z", dst_start_z);
        parameters1.put("nx", nx);
        parameters1.put("ny", ny);
        parameters1.put("nz", nz);
        long globalWorkSize1[] = {nx, ny, nz};
        clij2.execute(MyClij.class, "Kernels/getdata_periodic.cl", "copy_nxyz",
                globalWorkSize1, globalWorkSize1, parameters1);
    }

    public static void sum_up(CLIJ2 clij2, ClearCLImageInterface src, ClearCLImageInterface dst,
                              int nx, int ny, int nz) {
        HashMap<String, Object> parameters1 = new HashMap<>();
        parameters1.put("src", src);
        parameters1.put("dst", dst);
        long globalWorkSize1[] = {nx, ny, nz};
        clij2.execute(MyClij.class, "Kernels/sum_up.cl", "sum_up",
                globalWorkSize1, globalWorkSize1, parameters1);
    }


    public static void mul_add(CLIJ2 clij2, ClearCLImageInterface src,
                               int nx, int ny, int nz,
                               double alpha_cs, double delta) {
        HashMap<String, Object> parameters1 = new HashMap<>();
        parameters1.put("src", src);
        parameters1.put("mul", (float)alpha_cs);
        parameters1.put("add", (float)delta);

        long globalWorkSize1[] = {nx, ny, nz};
        clij2.execute(MyClij.class, "Kernels/mul_add.cl", "mul_add",
                globalWorkSize1, globalWorkSize1, parameters1);
    }






    public static void doTransform3D_clij(CLIJ2 clij2,
                                             ClearCLImageInterface image_stack_clij,
                                             ClearCLImageInterface Low_clij,
                                             ClearCLImageInterface buffer_clij,
                                             final int[] iterations, int nx, int ny, int nz) {
        assertDifferent(image_stack_clij, buffer_clij);
        int iterationx = iterations[0];
        int iterationy = iterations[1];
        int iterationz = iterations[2];
        int x_halfsize = nx;
        int y_halfsize = ny;
        int z_halfsize = nz;			//temporary variables used in the analysis
        final int inx=nx, iny=ny, inz=nz;
        int iteration = Math.max(iterationx, Math.max(iterationy, iterationz));
        for (int step = 1; step <= iteration; step++) {
            if (step <= iterationx) {
                x_halfsize = nx / 2;
                MyClij.haar_wavelet(clij2, image_stack_clij, buffer_clij, 0, nx, ny, nz, inx, iny, inz);
            }
            if (step <= iterationy) {
                y_halfsize = ny / 2;
                MyClij.haar_wavelet(clij2, image_stack_clij, buffer_clij, 1, nx, ny, nz, inx, iny, inz);
            }
            if (Low_clij != null) {
                HashMap<String, Object> parameters1 = new HashMap<>();
                parameters1.put("src", image_stack_clij);
                parameters1.put("dst", Low_clij);
                parameters1.put("src_start_x", 0);
                parameters1.put("src_start_y", 0);
                parameters1.put("src_start_z", 0);
                parameters1.put("dst_start_x", x_halfsize);
                parameters1.put("dst_start_y", y_halfsize);
                parameters1.put("dst_start_z", 0);
                long globalWorkSize1[] = {x_halfsize, y_halfsize, z_halfsize};
                clij2.execute(MyClij.class, "Kernels/getdata.cl", "copy_nxyz",
                        globalWorkSize1, globalWorkSize1, parameters1);
            }
            if (step <= iterationx) {
                nx = x_halfsize;
            }
            if (step <= iterationy) {
                ny = y_halfsize;
            }
            if (step <= iterationz) {
                nz = z_halfsize;
            }
        }
    }
    final static protected void doInverse3D_clij(CLIJ2 clij2, ClearCLImageInterface output_clij, ClearCLImageInterface buffer_clij,
                                                final int[] iterations, int nx, int ny, int nz) {

        int iterationx = iterations[0];
        int iterationy = iterations[1];
        int iterationz = iterations[2];
        int inx=nx, iny=ny, inz=nz;
        int x_size = nx / ((int) Math.pow((double) 2, (double) iterationx));
        int y_size = ny / ((int) Math.pow((double) 2, (double) iterationy));
        int z_size = nz / ((int) Math.pow((double) 2, (double) iterationz));

//        int x_halfsize = nx;
//        int y_halfsize = ny;
        int z_halfsize = nz;			//temporary variables used in the analysis
        int x_doublesize = (x_size == nx) ? nx : 2 * x_size;
        int y_doublesize = (y_size == ny) ? ny : 2 * y_size;
        int z_doublesize = (z_size == nz) ? nz : 2 * z_size;
        int iteration = Math.max(iterationx, Math.max(iterationy, iterationz));

        for (int step = iteration; step >= 1; step--) {
//            if (step <= iterationz) {
//                int kk=1;
//            }

            if (step <= iterationy) {
                MyClij.haar_wavelet(clij2, output_clij, buffer_clij, 3, x_doublesize, y_doublesize, z_doublesize, inx, iny, inz);
            }

            if (step <= iterationx) {
                MyClij.haar_wavelet(clij2, output_clij, buffer_clij, 2, x_doublesize, y_doublesize, z_doublesize, inx, iny, inz);
            }

            // reduce size for next iteration
            if (step <= iterationx) {
                x_size = x_doublesize;
                x_doublesize = (x_size == nx) ? nx : 2 * x_size;
            }
            if (step <= iterationy) {
                y_size = y_doublesize;
                y_doublesize = (y_size == ny) ? ny : 2 * y_size;
            }
//            if (step <= iterationz) {
//                nz = z_halfsize;
//            }

        } // go to next iteration

    }


    public static void show_clij_data(CLIJ2 clij2,
                                      ClearCLImageInterface image_stack_clij,
                                      String str
                                      ){

        ImagePlus F_imp = clij2.pull(image_stack_clij);
        ImageStack F_ist = F_imp.getStack();
        Object[] F_oj = F_ist.getImageArray();
        int inx = F_ist.getWidth();
        int iny = F_ist.getHeight();

        float[][] img_gpu = new float[inx][iny];
        float[] floats = (float[])F_oj[0];
        for (int i=0; i<inx; i++){
            for (int j=0; j<iny; j++){
//                img_gpu[i][j] = floats[i*iny+j];
                img_gpu[i][j] = floats[j*inx+i];
            }
        }
        ImageWare img_gpu_iw = Builder.create(inx, iny, 1, ImageWare.FLOAT);
        img_gpu_iw.putXY(0,0,0, img_gpu);
        img_gpu_iw.show(str);

    }

    public static void get_clij_data(CLIJ2 clij2,
                                      ClearCLImageInterface image_stack_clij,
                                      String str
    ){

        ImagePlus F_imp = clij2.pull(image_stack_clij);
        ImageStack F_ist = F_imp.getStack();
        Object[] F_oj = F_ist.getImageArray();
        int inx = F_ist.getWidth();
        int iny = F_ist.getHeight();

        float[][] img_gpu = new float[inx][iny];
        float[] floats = (float[])F_oj[0];
        for (int i=0; i<inx; i++){
            for (int j=0; j<iny; j++){
//                img_gpu[i][j] = floats[i*iny+j];
                img_gpu[i][j] = floats[j*inx+i];
            }
        }
        ImageWare img_gpu_iw = Builder.create(inx, iny, 1, ImageWare.FLOAT);
        img_gpu_iw.putXY(0,0,0, img_gpu);
        img_gpu_iw.show(str);

    }



    public static void haar_wavelet(CLIJ2 clij2, ClearCLImageInterface src, ClearCLImageInterface dst,
                                    int dir,
                                    int nx, int ny, int nz,
                                    int inx, int iny, int inz) {
        HashMap<String, Object> parameters = new HashMap<>();
        parameters.put("src", src);
        parameters.put("dst", dst);
        parameters.put("nx", nx);
        parameters.put("ny", ny);
        parameters.put("nz", nz);
        parameters.put("inx", inx);
        parameters.put("iny", iny);
        parameters.put("inz", inz);
        long globalWorkSize[] = {nx, ny, nz};
        String funcName = null;
        if (dir==0){
            funcName = "haar_wavelet_x";
        }else if (dir==1){
            funcName = "haar_wavelet_y";
        }else if (dir==2){
            funcName = "ihaar_wavelet_x";
        }else if (dir==3){
            funcName = "ihaar_wavelet_y";
        }
        clij2.execute(MyClij.class, "Kernels/wavelet.cl", funcName,
                globalWorkSize, globalWorkSize, parameters);

        HashMap<String, Object> parameters1 = new HashMap<>();
        parameters1.put("src", dst);
        parameters1.put("dst", src);
        parameters1.put("src_start_x", 0);
        parameters1.put("src_start_y", 0);
        parameters1.put("src_start_z", 0);
        parameters1.put("dst_start_x", 0);
        parameters1.put("dst_start_y", 0);
        parameters1.put("dst_start_z", 0);
        long globalWorkSize1[] = {nx, ny, nz};
        clij2.execute(MyClij.class, "Kernels/wavelet.cl", "copy_nxyz",
                globalWorkSize1, globalWorkSize1, parameters1);
    }

    public static void get_data(CLIJ2 clij2, ClearCLImageInterface src, ClearCLImageInterface dst,
                                int src_start_x, int src_start_y, int src_start_z,
                                int dst_start_x, int dst_start_y, int dst_start_z,
                                int nx, int ny, int nz) {
        HashMap<String, Object> parameters1 = new HashMap<>();
        parameters1.put("src", src);
        parameters1.put("dst", dst);
        parameters1.put("src_start_x", src_start_x);
        parameters1.put("src_start_y", src_start_y);
        parameters1.put("src_start_z", src_start_z);
        parameters1.put("dst_start_x", dst_start_x);
        parameters1.put("dst_start_y", dst_start_y);
        parameters1.put("dst_start_z", dst_start_z);
        long globalWorkSize1[] = {nx, ny, nz};
        clij2.execute(MyClij.class, "Kernels/getdata.cl", "copy_nxyz",
                globalWorkSize1, globalWorkSize1, parameters1);
    }



    final static protected double[] multiply(double[] array1, double[] array2) {
        double[] output = new double[array1.length];
        for (int i = 0; i < array1.length; i++) {
            output[i] = array1[i] * array2[i];
        }
        return output;
    }
    final static protected void multiply(double[] array1, double value, double[] output) {
        for (int i = 0; i < array1.length; i++) {
            output[i] = array1[i] * value;
        }
    }

    final static public void doMultiframePURELET_clij(CLIJ2 clij2, ClearCLBuffer image_stack_clij, ClearCLBuffer Low_clij, ClearCLBuffer output_clij,
                                                       int[] iterations, double[] sigma, int current,
                                                       int inx, int iny, int inz) {
//        MyClij.show_clij_data(clij2, Low_clij, "Low_clij");
        final boolean offsets[][] = {
                /* LLL */{false, false, false},
                /* LLH */ {false, false, true},
                /* LHL */ {false, true, false},
                /* LHH */ {false, true, true},
                /* HLL */ {true, false, false},
                /* HLH */ {true, false, true},
                /* HHL */ {true, true, false},
                /* HHH */ {true, true, true}
        };
        int nx = inx;
        int ny = iny;
        int nz = inz;
        int x_halfsize = inx;
        int y_halfsize = iny;
        int z_halfsize = inz;
        double[] sigma2 = multiply(sigma, sigma);
        double factor = 1.0;

//        double b1 = System.currentTimeMillis();
        int iters = Math.max(iterations[0], Math.max(iterations[1], iterations[2]));
//        int iters = 2;
        for (int i = 0; i < iters; i++) { // 每一个level
            factor = 1.0;
            if (iterations[0] > i) {
                x_halfsize = nx / 2; factor = 2.0 * factor;
            }
            if (iterations[1] > i) {
                y_halfsize = ny / 2; factor = 2.0 * factor;
            }
            if (iterations[2] > i) {
                z_halfsize = nz / 2; factor = 2.0 * factor;
            }
            multiply(sigma2, factor, sigma2); //得到当前level的sigma

            long[] dim_sub = {x_halfsize, y_halfsize, z_halfsize};
            ClearCLBuffer sub_low_clij = clij2.create(dim_sub);
            ClearCLBuffer sub_high_clij = clij2.create(dim_sub);
            long[] dim_sub1 = {x_halfsize, y_halfsize, 1};
            ClearCLBuffer sub_out_clij = clij2.create(dim_sub1);

            MyClij.get_data(clij2, Low_clij, sub_low_clij,
                    nx - x_halfsize, ny - y_halfsize, nz - z_halfsize,
                    0,0,0,
                    x_halfsize, y_halfsize, z_halfsize);
//            MyClij.show_clij_data(clij2, Low_clij, "Low_clij");

            int offset_x, offset_y;
            for (int j = 1; j < 8; j++) {
                if (!((iterations[0] <= i && offsets[j][0])  || (iterations[1] <= i && offsets[j][1])  || (iterations[2] <= i && offsets[j][2]))) {
                    offset_x = (offsets[j][0]) ? x_halfsize : 0 ;
                    offset_y = (offsets[j][1]) ? y_halfsize : 0 ;

                    MyClij.get_data(clij2, image_stack_clij, sub_high_clij,
                            offset_x, offset_y, 0,
                            0,0,0,
                            x_halfsize, y_halfsize, z_halfsize);

                    MyClij.PURELETmultishrink_clij(clij2, sub_low_clij, sub_high_clij, sub_out_clij,
                            offsets[j], sigma2, current,
                            x_halfsize, y_halfsize, z_halfsize);

                    MyClij.get_data(clij2, sub_out_clij, output_clij,
                            0,0,0,
                            offset_x, offset_y, 0,
                            x_halfsize, y_halfsize, 1);
                }
            }
            nx = x_halfsize;
            ny = y_halfsize;
            nz = z_halfsize;
            sub_low_clij.close();
            sub_high_clij.close();
            sub_out_clij.close();
        }

        long[] dim_sub1 = {x_halfsize, y_halfsize, z_halfsize};
        ClearCLBuffer low_clij = clij2.create(dim_sub1);
        long[] dim_sub2 = {x_halfsize, y_halfsize, 1};
        ClearCLBuffer out_clij = clij2.create(dim_sub2);

        MyClij.get_data(clij2, image_stack_clij, low_clij,
                0, 0, 0,
                0,0,0,
                x_halfsize, y_halfsize, z_halfsize);

        MyClij.DenoiseLowpass_clij(clij2, low_clij, out_clij,
                sigma2, current, x_halfsize, y_halfsize, z_halfsize);

        MyClij.get_data(clij2, out_clij, output_clij,
                0,0,0,
                0, 0, 0,
                x_halfsize, y_halfsize, 1);
        low_clij.close();
        out_clij.close();
    }










    public static void DenoiseLowpass_clij(CLIJ2 clij2, ClearCLBuffer low_clij, ClearCLBuffer out_clij,
                                           double[] sigma2, int currentFrame, int nx, int ny, int nz) {
        int N = nx * ny, K = nz;

        long[] dim_sub = {nz, nx, ny};
        ClearCLBuffer F_clij = clij2.create(dim_sub);
        ClearCLBuffer tmp_clij = clij2.create(dim_sub);


        /** get F and tmp */
        float[] sigma2_float = new float[sigma2.length];
        for (int k=0; k<sigma2.length; k++){
            sigma2_float[k] = (float)sigma2[k];
        }
        ClearCLBuffer sigma2_clij = clij2.pushArray(sigma2_float, sigma2_float.length, 1,1);
        HashMap<String, Object> parameters1 = new HashMap<>();
        parameters1.put("low", low_clij);
        parameters1.put("F", F_clij);
        parameters1.put("tmp", tmp_clij);
        parameters1.put("sigma2", sigma2_clij); //2
        parameters1.put("currentFrame", currentFrame);
        parameters1.put("nx", nx);
        parameters1.put("ny", ny);
        parameters1.put("nz", nz);
        long globalWorkSize1[] = {nx, ny};
        clij2.execute(MyClij.class, "Kernels/shrink_lowpass.cl", "shrink", globalWorkSize1, globalWorkSize1, parameters1);

//        Object F_oj = clij2.pullMatXYZ(F_clij);
//        Object tmp_oj = clij2.pullMatXYZ(tmp_clij);

        /** get C_float */
        int t = nx*ny;
        while (t > 1) {
            int m = t / 2;
            int n = (t + 1) / 2; // The fiddle with "m" and "n" is to handle non-power-of-two arrays.
            HashMap<String, Object> parameters5 = new HashMap<>();
            parameters5.put("tmp", tmp_clij);
            parameters5.put("K", K);
            parameters5.put("ny", ny);
            parameters5.put("offset", n);
            long globalWorkSize5[] = {K, m, 1};
            clij2.execute(MyClij.class, "Kernels/getC.cl", "getC", globalWorkSize5, globalWorkSize5, parameters5);
            t = n;
        }
        long[] dim_sub6 = {K, 1, 1};
        ClearCLBuffer C_clij = clij2.create(dim_sub6);
        MyClij.get_data(clij2, tmp_clij, C_clij,
                0,0,0,
                0, 0, 0,
                K, 1, 1);
        tmp_clij.close();
        ImagePlus C_imp = clij2.pull(C_clij);
        C_clij.close();
        ImageStack C_ist = C_imp.getStack();
        Object[] C_oj = C_ist.getImageArray();
        float[] C_float = (float[]) C_oj[0];


        /** get M_float */
        long[] dim_sub2 = {K*K, nx*ny, 1};
        ClearCLBuffer M_N = clij2.create(dim_sub2);
        HashMap<String, Object> parameters3 = new HashMap<>();
        parameters3.put("F", F_clij);
        parameters3.put("MN", M_N);
        parameters3.put("K", K);
        parameters3.put("nx", nx);
        parameters3.put("ny", ny);
        long globalWorkSize3[] = {K*K, nx, ny};
        clij2.execute(MyClij.class, "Kernels/getM.cl", "getMN", globalWorkSize3, globalWorkSize3, parameters3);
        t = nx*ny;
        while (t > 1) {
            int m = t / 2;
            int n = (t + 1) / 2; // The fiddle with "m" and "n" is to handle non-power-of-two arrays.
            HashMap<String, Object> parameters4 = new HashMap<>();
            parameters4.put("MN", M_N);
            parameters4.put("K", K);
            parameters4.put("offset", n);
            long globalWorkSize4[] = {K, K, m};
            clij2.execute(MyClij.class, "Kernels/fold.cl", "foldKernel", globalWorkSize4, globalWorkSize4, parameters4);
            t = n;
        }
        long[] dim_sub5 = {K*K, 1, 1};
        ClearCLBuffer M_clij = clij2.create(dim_sub5);
        MyClij.get_data(clij2, M_N, M_clij,
                0,0,0,
                0, 0, 0,
                K*K, 1, 1);
        M_N.close();
        ImagePlus M_imp = clij2.pull(M_clij);
        M_clij.close();
        ImageStack M_ist = M_imp.getStack();
        Object[] M_oj = M_ist.getImageArray();
        float[] M_float = (float[]) M_oj[0]; // have not divided by N



        Matrix C = new Matrix(nz, 1, 0.0);
        for (int c = 0; c < nz; c++) {
            C.set(c, 0, C_float[c]);
        }
        Matrix M = new Matrix(nz, nz, 0.0);
        for (int k1 = 0; k1 < nz; k1++) {
            for (int k2 = 0; k2 < nz; k2++) {
                M.set(k1, k2, M_float[k1*nz+k2] / N);
            }
        }
        //----------------------------------------------------------------------------------
        double lambda = (double) 1e-4;
        Matrix I = Matrix.identity(nz, nz);
        Matrix Mt = M.transpose();
        M = Mt.times(M).plus(I.times(lambda));
        //Matrix A = M.inverse().times(Mt).times(C);

        LUDecomposition LU = new LUDecomposition(M);
        Matrix A = Matrix.identity(K,K);
        if(!LU.isNonsingular()){
            IJ.log("LET parameters matrix is singular. Denoising results might be unreliable.");
            A = A.times(0);
        }
        else{
//            A = LU.solve(Mt.times(C));
            A = M.solve(Mt.times(C));
        }

        /** get LET */
        float[] A_float = new float[K];
        for (int i = 0; i < K; i++){
            A_float[i] = (float) A.get(i, 0);
        }
        ClearCLBuffer A_clij = clij2.pushArray(A_float, K, 1, 1);

        /** get output */
        HashMap<String, Object> parameters7 = new HashMap<>();
        parameters7.put("F", F_clij);
        parameters7.put("A", A_clij);
        parameters7.put("out", out_clij);
        parameters7.put("K", K);
        long globalWorkSize7[] = {nx, ny};
        clij2.execute(MyClij.class, "Kernels/getOutputLow.cl", "getOutput", globalWorkSize7, globalWorkSize7, parameters7);

//        Object out_oj = clij2.pullMatXYZ(out_clij);

        A_clij.close();
        F_clij.close();


    }

    final static public void PURELETmultishrink_clij(CLIJ2 clij2, ClearCLBuffer sub_low_clij, ClearCLBuffer sub_high_clij, ClearCLBuffer sub_out_clij,
        boolean[] offsets, double[] sigma2, int currentFrame, int nx, int ny, int nz) {
        int N = nx * ny, K = 6 * nz;

        long[] dim_sub = {nz*6, nx, ny};
        ClearCLBuffer F_clij = clij2.create(dim_sub);
        ClearCLBuffer tmp_clij = clij2.create(dim_sub);

        long[] dim_sub1 = {nx, ny, nz};
        ClearCLBuffer parent_clij = clij2.create(dim_sub1);

        /** get parent */
        int offset_arr[] = {(offsets[0]?1:0), (offsets[1]?1:0), (offsets[2]?1:0)};
        int offset_flag = (offsets[0]?1:0)*2 + (offsets[1]?1:0);
        HashMap<String, Object> parameters = new HashMap<>();
        parameters.put("low", sub_low_clij);
        parameters.put("dst", parent_clij);
        parameters.put("offset_flag", offset_flag); //1
        parameters.put("nx", nx);
        parameters.put("ny", ny);
        parameters.put("nz", nz);
        long globalWorkSize[] = {nx, ny, nz};
        clij2.execute(MyClij.class, "Kernels/shrink_highpass.cl", "get_parent", globalWorkSize, globalWorkSize, parameters);

        /** get F and tmp */
        float[] sigma2_float = new float[sigma2.length];
        for (int k=0; k<sigma2.length; k++){
            sigma2_float[k] = (float)sigma2[k];
        }
        ClearCLBuffer sigma2_clij = clij2.pushArray(sigma2_float, sigma2_float.length, 1,1);
        HashMap<String, Object> parameters1 = new HashMap<>();
        parameters1.put("high", sub_high_clij);
        parameters1.put("low", sub_low_clij);
        parameters1.put("parent", parent_clij);
        parameters1.put("F", F_clij);
        parameters1.put("tmp", tmp_clij);
        parameters1.put("sigma2", sigma2_clij); //2
        parameters1.put("currentFrame", currentFrame);
        parameters1.put("nx", nx);
        parameters1.put("ny", ny);
        parameters1.put("nz", nz);
        long globalWorkSize1[] = {nx, ny};
        clij2.execute(MyClij.class, "Kernels/shrink.cl", "shrink_high", globalWorkSize1, globalWorkSize1, parameters1);
        parent_clij.close();

        /** get C_float */
        int t = nx*ny;
        while (t > 1) {
            int m = t / 2;
            int n = (t + 1) / 2; // The fiddle with "m" and "n" is to handle non-power-of-two arrays.
            HashMap<String, Object> parameters5 = new HashMap<>();
            parameters5.put("tmp", tmp_clij);
            parameters5.put("K", K);
            parameters5.put("ny", ny);
            parameters5.put("offset", n);
            long globalWorkSize5[] = {K, m, 1};
            clij2.execute(MyClij.class, "Kernels/getC.cl", "getC", globalWorkSize5, globalWorkSize5, parameters5);
            t = n;
        }
        long[] dim_sub6 = {K, 1, 1};
        ClearCLBuffer C_clij = clij2.create(dim_sub6);
        MyClij.get_data(clij2, tmp_clij, C_clij,
                0,0,0,
                0, 0, 0,
                K, 1, 1);
        tmp_clij.close();
        ImagePlus C_imp = clij2.pull(C_clij);
        ImageStack C_ist = C_imp.getStack();
        Object[] C_oj = C_ist.getImageArray();
        float[] C_float = (float[]) C_oj[0];

        /** get M_float */
        long[] dim_sub2 = {K*K, nx*ny, 1};
        ClearCLBuffer M_N = clij2.create(dim_sub2);
        HashMap<String, Object> parameters3 = new HashMap<>();
        parameters3.put("F", F_clij);
        parameters3.put("MN", M_N);
        parameters3.put("K", K);
        parameters3.put("nx", nx);
        parameters3.put("ny", ny);
        long globalWorkSize3[] = {K*K, nx, ny};
        clij2.execute(MyClij.class, "Kernels/getM.cl", "getMN", globalWorkSize3, globalWorkSize3, parameters3);
        t = nx*ny;
        while (t > 1) {
            int m = t / 2;
            int n = (t + 1) / 2; // The fiddle with "m" and "n" is to handle non-power-of-two arrays.
            HashMap<String, Object> parameters4 = new HashMap<>();
            parameters4.put("MN", M_N);
            parameters4.put("K", K);
            parameters4.put("offset", n);
            long globalWorkSize4[] = {K, K, m};
            clij2.execute(MyClij.class, "Kernels/fold.cl", "foldKernel", globalWorkSize4, globalWorkSize4, parameters4);
            t = n;
        }
        long[] dim_sub5 = {K*K, 1, 1};
        ClearCLBuffer M_clij = clij2.create(dim_sub5);
        MyClij.get_data(clij2, M_N, M_clij,
                0,0,0,
                0, 0, 0,
                K*K, 1, 1);
        M_N.close();
        ImagePlus M_imp = clij2.pull(M_clij);
        ImageStack M_ist = M_imp.getStack();
        Object[] M_oj = M_ist.getImageArray();
        float[] M_float = (float[]) M_oj[0];

        /** get A */
        int k;
        double lambda = 1e-4;
        int[] indices = new int[K];
        Matrix M = new Matrix(K, K, 0.0);
        Matrix C = new Matrix(K, 1, 0.0);
        k = 0;
        for (int k1 = 0; k1 < K; k1++) {
            C.set(k1, 0, C_float[k1]);
            if (C_float[k1] >= 0) {
                indices[k] = k1;
                k++;
            }
            for (int k2 = 0; k2 < K; k2++) {
                M.set(k1, k2, M_float[k1*K+k2] / N);
            }
        }
        //----------------------------------------------------------------------------------
        K = k;
        Matrix Mr = new Matrix(K, K, 0.0);
        Matrix Cr = new Matrix(K, 1, 0.0);
        for (int k1 = 0; k1 < K; k1++) {
            Cr.set(k1, 0, C.get(indices[k1], 0));
            for (int k2 = 0; k2 < K; k2++) {
                Mr.set(k1, k2, M.get(indices[k1], indices[k2]));
            }
        }
        M = Mr.copy();
        C = Cr.copy();
        //----------------------------------------------------------------------------------
        Matrix I = Matrix.identity(K, K);
        Matrix Mt = M.transpose();
        Matrix MtC_ = Mt.times(C);
        M = Mt.times(M).plus(I.times(lambda));
        LUDecomposition LU = new LUDecomposition(M);
        Matrix A = Matrix.identity(K,K);
        if(!LU.isNonsingular()){
//            IJ.log("LET parameters matrix is singular. Denoising results might be unreliable.");
            A = A.times(0);
        }
        else{
            A = LU.solve(MtC_);
        }

        /** get LET */
        if (K != 0){
            float[] A_float = new float[K];
            float[] indices_float = new float[K];
            for (int i = 0; i < K; i++){
                A_float[i] = (float) A.get(i, 0);
                indices_float[i] = (float) indices[i];
            }
            ClearCLBuffer A_clij = clij2.pushArray(A_float, K, 1, 1);
            ClearCLBuffer indices_clij = clij2.pushArray(indices_float, K, 1, 1);

            /** get output */
            HashMap<String, Object> parameters7 = new HashMap<>();
            parameters7.put("F", F_clij);
            parameters7.put("A", A_clij);
            parameters7.put("out", sub_out_clij);
            parameters7.put("indices", indices_clij);
            parameters7.put("K", K);
            long globalWorkSize7[] = {nx, ny};
            clij2.execute(MyClij.class, "Kernels/getOutput.cl", "getOutput", globalWorkSize7, globalWorkSize7, parameters7);

//            Object F_oj = clij2.pullMatXYZ(F_clij);
//            Object A_oj = clij2.pullMatXYZ(A_clij);
//            Object sub_out_oj = clij2.pullMatXYZ(sub_out_clij);
//            Object indices_oj = clij2.pullMatXYZ(indices_clij);
            A_clij.close();
            indices_clij.close();
        }
        F_clij.close();
        tmp_clij.close();
    }
}
