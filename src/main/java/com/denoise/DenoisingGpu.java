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
package com.denoise;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.plugin.Memory;
import imageware.Builder;
import imageware.ImageWare;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij.coremem.enums.NativeTypeEnum;
import net.haesleinhuepf.clij2.CLIJ2;

import java.util.Random;
import java.util.Vector;

/*
 * Main class for multidimensional denoising of mixed Poisson-Gaussian noise.
 *
 * @author 	Swiss Federal Institute of Technology Lausanne
 *			Biomedical Imaging Group
 *
 * @version 1.0
 *
 */
public class DenoisingGpu {

    private ImageWare input, output;
    private double[] Alpha, Delta, Sigma;
    private int nx, ny, nz;
    private int CYCLESPIN, MULTIFRAME;
    private boolean FRAMEWISE = false;
    private boolean STOP = false;
    private boolean LOG = false;
    private boolean SUCCESS = false;
    private int MAX_THREAD = 20;
    private double MAX_MEMORY, IMAGE_MEMORY;
    private int[] CURRENT = new int[1];
    private Memory memo = new Memory();
    private double progress = 0;	// Status of the progress from 0 to 100.

    public CLIJ2 clij2;
    public ClearCLBuffer output_clij;
    /**
     * Constructor of the class Denoising.
     *
     * @param input         input data to be denoised (3D ImageWare object)
     * @param Alpha         double array containing the estimated detector gain 
     *                      for each frame/slice.
     * @param Delta         double array containing the estimated detector of-
     *                      fset for each frame/slice.
     * @param Sigma         double array containing the estimated AWGN standard
     *                      deviation for each frame/slice.
     * @param FRAMEWISE     true->Framewise noise parameters estimation
     *                      false->Global noise parameters estimation.
     * @param CYCLESPIN     number of cycle-spins (CS>0). A high value of CS
     *                      yields a high-quality denoising result, but the
     *                      computation time linearly increases with CS.
     * @param MULTIFRAME    number of adjacent frames/slices to be considered 
     *                      for multi-frame/slices denoising. MF>0 must be odd.
     *
     */
    public DenoisingGpu(ImageWare input,
            double[] Alpha,
            double[] Delta,
            double[] Sigma,
            boolean FRAMEWISE,
            int CYCLESPIN,
            int MULTIFRAME,
            CLIJ2 clij2) {
        this.input = input;
        this.Alpha = Alpha;
        this.Delta = Delta;
        this.Sigma = Sigma;
        this.FRAMEWISE = FRAMEWISE;
        this.CYCLESPIN = CYCLESPIN;
        this.MULTIFRAME = MULTIFRAME;
        this.clij2 = clij2;

        nx = input.getSizeX();
        ny = input.getSizeY();
        nz = input.getSizeZ();
        MAX_MEMORY = memo.maxMemory();//in bytes
//        IMAGE_MEMORY = nx * ny * nz * 4.0; //32 bit images in bytes
        IMAGE_MEMORY = nx * ny * nz * 4.0; //32 bit images in bytes
//        MAX_THREAD = (int) Math.max(Math.min(Math.round(0.25*MAX_MEMORY/IMAGE_MEMORY),nz),1);//Limit the number of concurrent threads
        MAX_THREAD = 1;
        System.out.printf("MAX_MEMORY: %f, IMAGE_MEMORY: %f, MAX_THREAD: %d \n", MAX_MEMORY, IMAGE_MEMORY, MAX_THREAD);
    }

    /**
     * Method for automatic (either global or framewise) noise parameters
     * estimation.
     */
    final public boolean estimateNoiseParameters() {
        ImageWare buffer = Builder.create(nx, ny, 1, ImageWare.DOUBLE);
        if (FRAMEWISE) {
            CURRENT[0] = 0;
            Vector<ParametersThread> V = new Vector<ParametersThread>();
            V.ensureCapacity(nz);
            ParametersThread pt = null;
            int k = 0, current_prev, k_prev;
            int Z = 0;
            while (Z < nz) {
                Z = Math.min(nz, Z + MAX_THREAD);
                k_prev = k;
                while (CURRENT[0] < Z) {
                    if (STOP) {
                        for (int j = k_prev; j < k; j++) {
                            pt = V.get(j);
                            pt.interrupt();
                        }
                        return false;
                    }
                    current_prev = CURRENT[0];
                    pt = new ParametersThread(input, Alpha, Delta, Sigma, CURRENT, SUCCESS, LOG);
                    pt.start();
                    if (LOG) {
                        IJ.log("Noise parameters estimation for frame no " + (current_prev + 1) + " lauched.");
                    }
                    pt.setPriority(Thread.MIN_PRIORITY);
                    V.insertElementAt(pt, k);
                    while (current_prev == CURRENT[0]) {
                        IJ.wait(1);
                    }
                    k++;
                }
                for (int i = k_prev; i < k; i++) {
                    pt = V.get(i);
                    if (STOP) {
                        for (int j = k_prev; j < k; j++) {
                            pt = V.get(j);
                            pt.interrupt();
                        }
                        return false;
                    }
                    try {
                        pt.join();
                        if (!pt.SUCCESS) {
                            for (int j = k_prev; j < k; j++) {
                                pt = V.get(j);
                                pt.interrupt();
                            }
                            return false;
                        }
                    } catch (Exception ex) {
                        IJ.log("Thread no " + i + ": " + pt.getState());
                        IJ.log("Error: thread cannot be ended.");
                    }
                    double p = 100 * (i + 1.0) / nz;
                    IJ.showStatus("Individual noise parameters estimation: " + IJ.d2s(p, 2) + "%");
                    System.out.println("Individual noise parameters estimation: " + IJ.d2s(p, 2) + "%");
                    progress = p;
                }
            }
            IJ.showStatus("");
            if (LOG) {
                IJ.log("---------------------------------------------------");
            }
        } else {
            IJ.showStatus("Global noise parameters estimation...");
            progress = 0;
            int c0 = (nz - 1) / 2;
            int zs = (int) Math.max(c0 - 1, 0);
            int ze = (int) Math.min(c0 + 1, nz - 1);
            buffer = Operations.averageSubStack(input, zs, ze);
            progress = 10;
            double[] noiseParams = Operations.estimateNoiseParams(buffer, 4);
            progress = 90;
            Alpha[0] = (ze - zs + 1) * noiseParams[0];
            Delta[0] = noiseParams[1];
            Sigma[0] = Math.sqrt(ze - zs + 1) * noiseParams[2];
            if (Alpha[0]<=0 || Sigma[0] < 0) {
                return false;
            }
            for (int i = 0; i < nz; i++) {
                Alpha[i] = Alpha[0];
                Delta[i] = Delta[0];
                Sigma[i] = Sigma[0];
            }
            progress = 100;
            IJ.showStatus("");
        }

        return true;
    }

    /**
     * Method for parallel denoising of multidimensional images.
     */


    final public void perform() {
        progress = 0;
        int[] iters = Operations.numOfDyadicIters(nx, ny);

//        CLIJ2 clij2 = CLIJ2.getInstance("AMD");
        CLIJ2 clij2 = this.clij2;
        String s = clij2.getGPUName();
        System.out.println(s);
//        IJ.log("Device: " + s);

        // cpu -> gpu
        ImagePlus input_imp = new ImagePlus("source", input.buildImageStack());
        ClearCLBuffer input_clij = clij2.push(input_imp);

        if (input.getDimension() == 2) {
            ClearCLBuffer input_clij_2d = clij2.push(input_imp);
            long[] dim_a = {nx, ny, 1};
            input_clij = clij2.create(dim_a, NativeTypeEnum.Short);
            MyClij.extend_dim(clij2, input_clij_2d, input_clij,
                    0, 0, 0,
                    0, 0, 0,
                    nx, ny, 1);
        }
        ClearCLBuffer output_clij = clij2.create(input_clij);

//        MyClij.show_clij_data(clij2, input_clij, "input_clij");

        // do processing in gpu
        long[] dim_oneframe = {nx, ny, 1};

        double t1 = System.currentTimeMillis();

        for (int current_p = 0; current_p < nz; current_p++) {
//        for (int current_p = 0; current_p < 10; current_p++) {
//            ClearCLBuffer out_sum_clij = clij2.create(dim_oneframe);
            float[][][] float_arr = new float[nx][ny][1];
            ClearCLBuffer out_sum_clij = clij2.pushMatXYZ(float_arr);


            int[] index = Operations.getAdjacentIndex(current_p, nz, MULTIFRAME);
            double[] bufferSigma = Operations.restrictArray(Operations.divide(Sigma, Alpha), index[0], index[1]);
            int len_stack = index[1] - index[0] + 1;


            long[] dim_sub = {nx, ny, len_stack};
            ClearCLBuffer coef_clij = clij2.create(dim_sub);
            ClearCLBuffer image_stack_clij = clij2.create(dim_sub);

            MyClij.get_data(clij2, input_clij, coef_clij,
                    0, 0, index[0],
                    0, 0, 0,
                    nx, ny, len_stack);

            MyClij.sub_div(clij2, coef_clij, nx, ny, len_stack, Alpha, Delta, index);
            MyClij.show_clij_data(clij2, coef_clij, "coef_clij");


            Random rand = new Random(0);
            for (int cs = 0; cs < CYCLESPIN; cs++) {
                int[] offset = new int[2];
                offset[0] = (int) (cs > 0 ? Math.floor(rand.nextDouble() * nx) : 0);
                offset[1] = (int) (cs > 0 ? Math.floor(rand.nextDouble() * ny) : 0);

                MyClij.get_data_periodic(clij2, coef_clij, image_stack_clij,
                        offset[0], offset[1], 0,
                        0, 0, 0,
                        nx, ny, len_stack);
//                MyClij.show_clij_data(clij2, image_stack_clij, "image_stack_clij");
                ClearCLBuffer out_clij = clij2.create(dim_oneframe);
                ClearCLBuffer out_shift_clij = clij2.create(dim_oneframe);
                ClearCLBuffer Low_clij = clij2.create(image_stack_clij);
                ClearCLBuffer tmp_clij = clij2.create(image_stack_clij);
                MyClij.doTransform3D_clij(clij2, image_stack_clij, Low_clij, tmp_clij, iters, nx, ny, len_stack);
                MyClij.doMultiframePURELET_clij(clij2, image_stack_clij, Low_clij, out_clij, iters, bufferSigma, index[2], nx, ny, len_stack);
                MyClij.doInverse3D_clij(clij2, out_clij, tmp_clij, iters, nx, ny, 1);
//                clij2.copy(image_stack_clij, out_clij);

                MyClij.get_data_periodic(clij2, out_clij, out_shift_clij,
                        nx - offset[0], ny - offset[1], 0,
                        0, 0, 0,
                        nx, ny, len_stack);
                MyClij.sum_up(clij2, out_shift_clij, out_sum_clij, nx, ny, 1);
                out_clij.close();
                out_shift_clij.close();
                Low_clij.close();
                tmp_clij.close();
            }
            MyClij.mul_add(clij2, out_sum_clij, nx, ny, 1, Alpha[current_p] / CYCLESPIN, Delta[current_p]);
//            MyClij.show_clij_data(clij2, out_sum_clij, "mul_add");


            MyClij.get_data(clij2, out_sum_clij, output_clij,
                    0, 0, 0,
                    0, 0, index[0] + index[2],
                    nx, ny, 1);

            coef_clij.close();
            image_stack_clij.close();
            out_sum_clij.close();

            double t2 = System.currentTimeMillis();
            System.out.printf("Denoising of frame no %d launched. \n", (current_p + 1));
            IJ.log("Denoising of frame no " + String.valueOf((current_p + 1)) + " done. ( " + String.valueOf((t2-t1)/1000) + "s )" );
            t1 = t2;
        }


//        MyClij.show_clij_data(clij2, output_clij, "output_clij");
        ImagePlus F_imp = clij2.pull(output_clij);
        ImageStack F_ist = F_imp.getStack();
        Object[] F_oj = F_ist.getImageArray();
        int inx = F_ist.getWidth();
        int iny = F_ist.getHeight();
        int inz = F_ist.getSize();
        input_clij.close();
        output_clij.close();

        output = Builder.create(nx, ny, nz, ImageWare.SHORT);
        short[][] img_gpu = new short[inx][iny];
        for (int k=0; k<nz; k++) {
            short[] floats = (short[]) F_oj[k];
            for (int i = 0; i < inx; i++) {
                for (int j = 0; j < iny; j++) {
                    img_gpu[i][j] = floats[j * inx + i];
                }
            }
            output.putXY(0,0, k, img_gpu);
        }

    }


    /**
     * Method to set the type of noise parmeters estimation (global or framewise).
     *
     * @param framewise  type of noise parmeters estimation:
     *                   true->Framewise.
     *                   false->Global.
     */
    final public void setFramewise(boolean framewise) {
        FRAMEWISE = framewise;
    }

    /**
     * Method to set the number of cycle-spin(s) for reducing potential blocking
     * artifacts.
     *
     * @param cyclespin  number of cycle-spins (cyclespin>0).
     */
    final public void setCycleSpins(int cyclespin) {
        CYCLESPIN = cyclespin;
    }

    /**
     * Method to interrupt the current denoising task.
     *
     * @param stop  true->Stop the current denoising task.
     */
    final public void setStop(boolean stop) {
        STOP = stop;
    }

    /**
     * Method to display some messages related to the current denoising task.
     *
     * @param log    true->Enable display.
     *               false->Disable display.
     */
    final public void setLog(boolean log) {
        LOG = log;
    }

    /**
     * Method to set the number of adjacent frames/slices for multi-frame/slices
     * denoising.
     *
     * @param multiframe  number of adjacent frames/slices for multi-frame/slices
     *                    denoising.
     */
    final public void setMultiFrame(int multiframe) {
        MULTIFRAME = multiframe;
    }

    /**
     * Method to set the value of the detector gain for each frame/slice.
     *
     * @param alpha  array containing the value of the detector gain for each frame/slice.
     */
    final public void setAlpha(double[] alpha) {
        for (int i = 0; i < nz; i++) {
            Alpha[i] = alpha[i];
        }
    }

    /**
     * Method to set the value of the detector offset for each frame/slice.
     *
     * @param delta  array containing the value of the detector offset for each frame/slice.
     */
    final public void setDelta(double[] delta) {
        for (int i = 0; i < nz; i++) {
            Delta[i] = delta[i];
        }
    }

    /**
     * Method to set the standard deviation of the additive-white-Gaussian-noise
     * (AWGN) for each frame/slice.
     *
     * @param sigma  array containing the standard deviation of the AWGN for each
     *               frame/slice.
     */
    final public void setSigma(double[] sigma) {
        for (int i = 0; i < nz; i++) {
            Sigma[i] = sigma[i];
        }
    }

    /**
     * Method to get the value of the detector gain for each frame/slice.
     *
     * @return  Alpha  array containing the value of the detector gain for each frame/slice.
     */
    final public double[] getAlpha() {
        return Alpha;
    }

    /**
     * Method to get the value of the detector offset for each frame/slice.
     *
     * @return Delta  array containing the value of the detector offset for each frame/slice.
     */
    final public double[] getDelta() {
        return Delta;
    }

    /**
     * Method to get the standard deviation of the additive-white-Gaussian-noise
     * (AWGN) for each frame/slice.
     *
     * @return Sigma  array containing the standard deviation of the AWGN for each
     *                frame/slice.
     */
    final public double[] getSigma() {
        return Sigma;
    }

    /**
     * Method to get the output of the denoising task.
     *
     * @return output  3D ImageWare containing the denoised image.
     */
    final public ImageWare getOutput() {
        return output;
    }

    /**
     * Method to get the maximum number of parallel threads launched.
     *
     * @return MAX_THREAD  the maximum number of parallel threads launched.
     */
    final public int getMaxThread() {
        return MAX_THREAD;
    }

    /**
     * Method to get the level of progression during the denoising process.
     *
     * @return progress level between 0 and 100.
     */
    final public double getProgress() {
        return progress;
    }

    /**
     * Method to detect whether the denoising process has been stopped by the user.
     *
     * @return STOP boolean.
     */
    final public boolean getStop() {
        return STOP;
    }
}