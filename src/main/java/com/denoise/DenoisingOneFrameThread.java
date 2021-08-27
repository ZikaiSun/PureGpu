package com.denoise;

//====================================================================
//
// Project: PureDenoise
// Package: denoise
// Class  : DenoisingThread
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

import ij.*;
import imageware.Builder;
import imageware.ImageWare;


import ij.IJ;
import ij.plugin.Duplicator;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;
import net.haesleinhuepf.clij2.plugins.Flip2D;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.integer.UnsignedShortType;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import net.haesleinhuepf.clij2.AbstractCLIJ2Plugin;
import net.haesleinhuepf.clij2.plugins.Flip2D;
import org.jocl.Pointer;
import org.jocl.Sizeof;
import org.jocl.cl_mem;
import static org.jocl.CL.*;


public class DenoisingOneFrameThread extends Thread  {

    ImageWare Coef, output;
    double Alpha, Delta;
    double[] Sigma;
    int[] current, iters;
    int pos;
    int filter = Filters.UNHAAR;
    double degree = 1.0, shift = 0.0;
    boolean LOG;
    int CYCLESPIN;

    protected DenoisingOneFrameThread(ImageWare Coef,
                                      ImageWare output,
                                      double Alpha,
                                      double Delta,
                                      double[] Sigma,
                                      int[] iters,
                                      int pos,
                                      int[] current,
                                      boolean LOG,
                                      int CYCLESPIN) {
        this.Coef = Coef;
        this.output = output;
        this.Alpha = Alpha;
        this.Delta = Delta;
        this.Sigma = Sigma;
        this.iters = iters;
        this.pos = pos;
        this.current = current;
        this.LOG = LOG;
        this.CYCLESPIN = CYCLESPIN;
    }

    @Override
    public void run() {
        int c = current[0];
        current[0] = current[0] + 1;
        int nx = Coef.getSizeX(), ny = Coef.getSizeY(), nz = Coef.getSizeZ();
        java.util.Random rand = new java.util.Random(0);

        ImageWare out = Builder.create(nx, ny, nz, ImageWare.DOUBLE);

        /** --CPU version-- */
        out = Builder.create(nx, ny, nz, ImageWare.DOUBLE);
        for (int cs = 0; cs < CYCLESPIN; cs++) {

            int[] offset = new int[2];
            offset[0] = (int) (cs > 0 ? Math.floor(rand.nextDouble() * nx) : 0);
            offset[1] = (int) (cs > 0 ? Math.floor(rand.nextDouble() * ny) : 0);

            double[][][] array = new double[nx][ny][Sigma.length];
            Coef.getBlockXYZ(offset[0], offset[1], 0, array, ImageWare.PERIODIC);

            double c1 = System.currentTimeMillis();

            ImageWare image_stack = Builder.create(array, ImageWare.DOUBLE);
            ImageWare Low = Builder.create(nx, ny, nz, ImageWare.DOUBLE);
            Operations.doTransform3D(image_stack, iters, filter, degree, shift, 0, Low);
            ImageWare buffer = Operations.doMultiframePURELET(image_stack, Low, iters, Sigma, pos);
            Operations.doInverse3D(buffer, iters, filter, degree, shift);

            double c2 = System.currentTimeMillis();
            System.out.println(c2-c1 + " ms (cpu)");

            double[][][] output_cpu = new double[nx][ny][1];
            buffer.getBlockXYZ(-offset[0], -offset[1], 0, output_cpu, ImageWare.PERIODIC);

            buffer.multiply(Alpha);
            buffer.add(Delta);
            array = new double[nx][ny][nz];
            buffer.getBlockXYZ(-offset[0], -offset[1], 0, array, ImageWare.PERIODIC);
            out.add(Builder.create(array, ImageWare.DOUBLE));
        }
        out.divide(CYCLESPIN);
//        out.show("out");

        output.putXY(0, 0, c, out);

        if (LOG) {
            IJ.log("Frame no " + (c + 1) + " denoised.");
        }

        output.show("output at frame");
    }
}
