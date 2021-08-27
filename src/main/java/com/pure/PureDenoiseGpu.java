package com.pure;

import com.denoisegui.DenoiseDialogGpu;

import ij.ImagePlus;
import ij.LookUpTable;
import ij.WindowManager;
import imageware.Builder;
import imageware.ImageWare;
import net.haesleinhuepf.clij.macro.CLIJMacroPlugin;
import net.haesleinhuepf.clij.macro.CLIJOpenCLProcessor;
import net.haesleinhuepf.clij.macro.documentation.OffersDocumentation;
import net.haesleinhuepf.clij2.AbstractCLIJ2Plugin;
import net.haesleinhuepf.clij2.CLIJ2;
import net.haesleinhuepf.clij2.utilities.HasAuthor;
import net.haesleinhuepf.clij2.utilities.HasClassifiedInputOutput;
import net.haesleinhuepf.clij2.utilities.HasLicense;
import net.haesleinhuepf.clij2.utilities.IsCategorized;
import org.scijava.plugin.Plugin;

import com.denoise.DenoisingGpu;
import com.denoise.Operations;

import ij.IJ;

import javax.swing.*;


/**
 * code modified from http://bigwww.epfl.ch/algorithms/denoise/
 * Author: Zikai Sun
 *         Apr 2021
 */
@Plugin(type = CLIJMacroPlugin.class, name = "CLIJ2_PureDenoiseGpu")
public class PureDenoiseGpu extends AbstractCLIJ2Plugin implements CLIJMacroPlugin, CLIJOpenCLProcessor, HasClassifiedInputOutput, IsCategorized {

    @Override

    public String getParameterHelpText() {
//        return "Image source, ByRef Image destination";
        return "";
    }

    @Override
    public boolean executeCL() {

//        CLIJ2 clij2 = CLIJ2.getInstance("AMD");
        CLIJ2 clij2 = getCLIJ2();
        String s = clij2.getGPUName();
//        IJ.log(s);

        new DenoiseDialogGpu(clij2);

        return true;
    }

    @Override
    public String getCategories() {
        return "Filter";
    }

    @Override
    public String getInputType() {
        return "Image";
    }

    @Override
    public String getOutputType() {
        return "Image";
    }

    public static int get_num(int a){
        int n = a -1;
        n |= n >>> 1;
        n |= n >>> 2;
        n |= n >>> 4;
        n |= n >>> 8;
        n |= n >>> 16;
        return n+1;
    }

    public static void main(final String... args)  {

        CLIJ2 clij2 = CLIJ2.getInstance("AMD");

        String titleIn = "";
        String titleOut = "";
        byte[] red, green, blue;
        int nx = 0;
        int ny = 0;
        int nz = 0;
        double Imax = 0.0, Imin = 0.0;
        ImageWare original = null;		// Input image
        ImageWare output = null;		// Output image
        boolean COLOR = false;
        int Nmin = 16;
        int nxe = 0;
        int nye = 0;
        int[] Ext = new int[2];
        double[] AlphaHat, DeltaHat, SigmaHat;	// estimated gain & offset of the detectors, AGWN std
        DenoisingGpu denoising;
        boolean FRAMEWISE = false;
        int CS = 1;
        int NBFRAME = 0;
        boolean LOG = false;
        double Alpha = 1.0, Delta = 0.0, Sigma = 0.0;

        IJ.open("/Users/jake3un/Documents/Noisy_Cell_Colony.jpg");

        double time = System.currentTimeMillis();
        ImagePlus impSource = WindowManager.getCurrentImage();
        int type = impSource.getType();
        titleIn = impSource.getTitle();
        titleOut = "Denoised-" + titleIn;
        LookUpTable lut = impSource.createLut();
        red = lut.getReds();
        green = lut.getGreens();
        blue = lut.getBlues();

        nx = impSource.getWidth();
        ny = impSource.getHeight();
        nz = impSource.getStackSize();

        Imax = impSource.getDisplayRangeMax();
        Imin = impSource.getDisplayRangeMin();

        if (type == ImagePlus.COLOR_RGB) {
            if (nz > 1) {
                IJ.showMessage("Note", "This version of the plugin does not handle color stacks.");
                return;
            } else {
                COLOR = true;
                nz = 3;
                original = Builder.create(nx, ny, nz, ImageWare.DOUBLE);
                ImageWare[] temp = Builder.createColors(impSource);
                original.putXY(0, 0, 0, temp[0]);
                original.putXY(0, 0, 1, temp[1]);
                original.putXY(0, 0, 2, temp[2]);
            }
        } else {
            COLOR = false;
            original = Builder.create(impSource);
        }
        if (nx < Nmin || ny < Nmin) {
            IJ.log("The size of your data is inapropriate.");
            return;
        }
        if (original == null) {
            IJ.log("Unable to create the data set.");
            return;
        }
        IJ.showStatus("Denoising in progress...");

        impSource = null;
        original = original.convert(ImageWare.SHORT);
        nz = original.getSizeZ();


        nxe = get_num(nx);
        nye = get_num(ny);

//        nxe = (int) (Math.ceil((double) nx / Nmin) * Nmin);
//        nye = (int) (Math.ceil((double) ny / Nmin) * Nmin);
        if (nxe != nx || nye != ny) {
            original = Operations.symextend2D(original, nxe, nye, Ext);
        } else {
            Ext[0] = 0;
            Ext[1] = 0;
        }
        AlphaHat = new double[nz];
        DeltaHat = new double[nz];
        SigmaHat = new double[nz];

        JRadioButton checkAutomatic;
        denoising  = new DenoisingGpu(original,AlphaHat,DeltaHat,SigmaHat, FRAMEWISE, CS, NBFRAME, clij2);

        denoising.setLog(LOG);
//        String argsEstimation[] = {"Auto", "Global"};
        String argsEstimation[] = {"Manual", "1", "178.673", "26.641"};

        for (int u = 0; u < argsEstimation.length; u++) {
            IJ.log("argsEstimation["+u+"]: "+argsEstimation[u]);
        }
        if (argsEstimation.length != 2 && argsEstimation.length != 4) {
            IJ.error("The estimation parameters is incorrect. Correct example:\"estimation=Auto Global\" or \"estimation=Manual 30.0 3.0 40.0\" ");
            return;
        }

        if (argsEstimation[0].toLowerCase().equals("auto")) {
            denoising.setFramewise(argsEstimation[1].toLowerCase().equals("individual"));
            denoising.estimateNoiseParameters();
        } else if (argsEstimation[0].toLowerCase().equals("manual")) {
            Alpha = (new Double(argsEstimation[1])).doubleValue();
            Delta = (new Double(argsEstimation[2])).doubleValue();
            Sigma = (new Double(argsEstimation[3])).doubleValue();
            for (int i = 0; i < nz; i++) {
                AlphaHat[i] = Alpha;
                DeltaHat[i] = Delta;
                SigmaHat[i] = Sigma;
            }
        } else {
            IJ.error("The estimation parameters is incorrect. Correct example:\"estimation=Auto Global\" or \"estimation=Manual 30.0 3.0 40.0\" ");
            return;
        }

        String argsParameters[] = {"3", "4"};
        for (int u = 0; u < argsParameters.length; u++) {
            IJ.log("argsParameters["+u+"]: "+argsParameters[u]);
        }
        if (argsParameters.length != 2) {
            IJ.error("The parameters is incorrect. Correct example:\"parameters=3 4\" ");
            return;
        }
        NBFRAME = Math.max(1, Math.min(nz, (int) (new Double(argsParameters[0])).doubleValue()));
        CS = Math.max(1, Math.min(10, (int) (new Double(argsParameters[1])).doubleValue()));
        denoising.setCycleSpins(CS);
        denoising.setMultiFrame(NBFRAME);

        denoising.perform();
        output = denoising.getOutput();
        if (nxe != nx || nye != ny) {
            output = Operations.crop2D(output, nx, ny, Ext);
        }
        output.show("Denoised " + titleIn);
        //----------------------------------------------------------------------------------
        IJ.log("-------------- SUMMARY --------------");
        IJ.log("Noise parameters used for denoising: \"" + titleIn + "\"");
        for (int i = 0; i < nz; i++) {
            IJ.log("Frame " + (i + 1) + ": Alpha = " + IJ.d2s(AlphaHat[i], 3) + " Delta = " + IJ.d2s(DeltaHat[i], 3) + " Sigma = " + IJ.d2s(SigmaHat[i], 3));
        }
        IJ.log("Number of adjacent frames = " + NBFRAME);
        IJ.log("Number of cycle-spin(s) = " + CS);
        IJ.log("Maximum number of concurrent threads: " + denoising.getMaxThread());
        IJ.log("The whole processing required " + IJ.d2s((System.currentTimeMillis() - time) / 1000.0, 2) + " s.");
        IJ.log("-------------------------------------");





//        PluginTemplate pd = new PluginTemplate();
//        pd.executeCL();
//        new DenoiseDialogGpu();

    }

}