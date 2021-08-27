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

import ij.IJ;
import imageware.Builder;
import imageware.ImageWare;


public class DenoisingThread extends Thread {

    ImageWare Coef, output;
    double Alpha, Delta;
    double[] Sigma;
    int[] current, iters;
    int pos;
    int filter = Filters.UNHAAR;
    double degree = 1.0, shift = 0.0;
    boolean LOG;

    protected DenoisingThread(ImageWare Coef,
            ImageWare output,
            double Alpha,
            double Delta,
            double[] Sigma,
            int[] iters,
            int pos,
            int[] current,
            boolean LOG) {
        this.Coef = Coef;
        this.output = output;
        this.Alpha = Alpha;
        this.Delta = Delta;
        this.Sigma = Sigma;
        this.iters = iters;
        this.pos = pos;
        this.current = current;
        this.LOG = LOG;
    }

    @Override
    public void run() {
        int c = current[0];
        current[0] = current[0] + 1;
        double c0 = System.currentTimeMillis();
        int nx = Coef.getSizeX(), ny = Coef.getSizeY(), nz = Coef.getSizeZ();
        ImageWare Low = Builder.create(nx, ny, nz, ImageWare.DOUBLE);

//        double c1 = System.currentTimeMillis();
        Operations.doTransform3D(Coef, iters, filter, degree, shift, 0, Low);
//        double c2 = System.currentTimeMillis();
        ImageWare buffer = Operations.doMultiframePURELET(Coef, Low, iters, Sigma, pos);
//        double c3 = System.currentTimeMillis();
        Operations.doInverse3D(buffer, iters, filter, degree, shift);
//        double c4 = System.currentTimeMillis();

        buffer.multiply(Alpha);
        buffer.add(Delta);
        output.putXY(0, 0, c, buffer);
        if (LOG) {
            IJ.log("Frame no " + (c + 1) + " denoised.");
        }

    }
}
