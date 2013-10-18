package Mycelium;

import AnaMorf.BatchAnalyser;
import IAClasses.FractalEstimator;
import IAClasses.ProgressDialog;
import IAClasses.Utils;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.Roi;
import ij.plugin.filter.GaussianBlur;
import ij.process.ByteProcessor;
import ij.process.FloatBlitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;

public class MyceliumGrower {

    public double radius = 25.0;
    private FloatProcessor densityField, nutrientField;
    public double gradSens = 5000.0, noise = 10.0, growThresh = 150.0;
    public ArrayList dataPoints;
    private static int popSize = 1;
    public File imageFolder;
    public ImageProcessor[] lacSurf, dsSurf, areaSurf;
    public double areaCurve[];
    public double ds_max = -Double.MAX_VALUE, lac_max = -Double.MAX_VALUE,
            area_max = -Double.MAX_VALUE;
    private Random R = new Random();

    public static void main(String args[]) {
        MyceliumGrower grower = new MyceliumGrower();
        float steps = 40000.0f;
//        Random rand = new Random();
//        double maxinc = 5.0;
//        grower.lacSurf = new ImageProcessor[10];
//        grower.dsSurf = new ImageProcessor[10];
//        grower.areaSurf = new ImageProcessor[10];
//        for (int i = 0; i < 10; i++) {
//            grower.lacSurf[i] = (new ImagePlus("C:\\Users\\Dave\\lac" + (i + 1) + ".tif")).getProcessor();
//            grower.dsSurf[i] = (new ImagePlus("C:\\Users\\Dave\\ds" + (i + 1) + ".tif")).getProcessor();
//            grower.areaSurf[i] = (new ImagePlus("C:\\Users\\Dave\\area" + (i + 1) + ".tif")).getProcessor();
//            grower.areaSurf[i].log();
//            ImageStatistics lacstats = grower.lacSurf[i].getStatistics();
//            ImageStatistics fracstats = grower.dsSurf[i].getStatistics();
//            ImageStatistics areastats = grower.areaSurf[i].getStatistics();
//            if (lacstats.max > grower.lac_max) {
//                grower.lac_max = lacstats.max;
//            }
//            if (fracstats.max > grower.ds_max) {
//                grower.ds_max = fracstats.max;
//            }
//            if (areastats.max > grower.area_max) {
//                grower.area_max = areastats.max;
//            }
//        }
//        grower.areaCurve = new double[grower.areaSurf[0].getHeight()];
//        Arrays.fill(grower.areaCurve, 0.0);
//        for (int i = 0; i < 10; i++) {
//            grower.lacSurf[i].multiply(1.0 / grower.lac_max);
//            grower.dsSurf[i].multiply(1.0 / grower.ds_max);
//            grower.areaSurf[i].multiply(1.0 / grower.area_max);
//            for (int y = 0; y < grower.areaSurf[i].getHeight(); y++) {
//                double sum = 0.0;
//                for (int x = 0; x < grower.areaSurf[i].getWidth(); x++) {
//                    sum += grower.areaSurf[i].getPixelValue(x, y);
//                }
//                grower.areaCurve[y] += sum / (grower.areaSurf[i].getWidth() * 10.0);
//            }
//        }
//        for (double h = 10+ maxinc * rand.nextDouble(); h < 100.0; h += maxinc * rand.nextDouble()) {
        for (int h = 50; h <= 50; h += 5) {
            grower.dataPoints = new ArrayList();
            for (int i = 0; i < popSize; i++) {
                ByteProcessor bp = new ByteProcessor(1200, 1200);
                bp.setColor(Color.white);
                bp.fill();
                grower.run(bp, i, h, steps);
            }
            File allResults = null;
            PrintWriter outputStreamAll = null;
            try {
                allResults = new File(grower.imageFolder + "\\..\\Results" + h
                        + "_" + steps + "_" + popSize + ".csv");
            } catch (Exception e) {
                IJ.error(e.toString());
            }
            try {
                outputStreamAll = new PrintWriter(new FileOutputStream(allResults));
            } catch (FileNotFoundException e) {
                IJ.error("Could not write to results file.");
            }
            outputStreamAll.println("Growth Unit:," + h + ",Filter Radius:," + grower.radius
                    + ",Gradient Sensitivity:," + grower.gradSens + ",Noise:," + grower.noise + "\n");
            outputStreamAll.println("Iteration,Number of Tips,,Total Length,,Perimeter Length,,Circularity,,Area,,Dbm,,Dbs,,Ds,,Dss,,Lacunarity,,Ds R^2,,Dss R^2,,Dst,,Dst R^2,,HGU Estimate,,N");
            for (int k = 0; k < grower.dataPoints.size(); k++) {
                double currentRes[][] = (double[][]) grower.dataPoints.get(k);
                outputStreamAll.print(k);
                for (int l = 0; l < 15; l++) {
                    double[] thisresult = getMeanAndSD(currentRes[l]);
                    outputStreamAll.print("," + thisresult[0] + "," + thisresult[1]);
                }
                outputStreamAll.print("," + currentRes[15][0] + "\n");
            }
            outputStreamAll.close();
        }
        System.exit(0);
    }

    /*
     * private boolean showDialog() { GenericDialog gd = new
     * GenericDialog(this.getClass().getName()); gd.addNumericField("HGU:", hgu,
     * 0, 6, "pixels"); gd.addNumericField("Iterations:", iterations, 0, 6, "
     * "); gd.showDialog(); if (gd.wasCanceled()) { return false; } hgu =
     * gd.getNextNumber(); iterations = gd.getNextNumber(); return true; }
     */
    public MyceliumGrower() {
    }

    public void run(ByteProcessor ip, int thisIter, double hgu, float maxLength) {
        int w = ip.getWidth(), h = ip.getHeight(), i, j, h0 = 0, h1, x0 = w / 2,
                y0 = h / 2, x, y, xlow = w / 2 - 1, xhigh = w / 2 + 1,
                ylow = h / 2 - 1, yhigh = h / 2 + 1;
        double apex;
        FractalEstimator boxCounter = new FractalEstimator();
        GaussianBlur blurrer = new GaussianBlur();
        DecimalFormat numFormat = new DecimalFormat("000");
        DecimalFormat decFormat = new DecimalFormat("0.000");
        int totalLength = 0;

        ip.setValue(0);
        ArrayList hyphae = new ArrayList();
        double angle = R.nextDouble() * 360;

        hyphae.add(new Hypha(x0, y0, angle, ip, hgu));
        h0++;
        /*
         * double angle2 = angle + 180.0; if (angle2 > 360) { angle2 -= 360.0; }
         * hyphae.add(new Hypha(x0, y0, angle2, ip, hgu)); h0++;
         */

        imageFolder = new File("c:\\Mycelium\\" + decFormat.format(hgu)
                + "_" + numFormat.format(maxLength));
        if (!imageFolder.exists()) {
            if (!imageFolder.mkdir()) {
                IJ.error("Failed to create image directory.");
            }
        }
        File results = null;
        PrintWriter outputStream = null;
        try {
            results = new File(imageFolder + "\\Results_" + thisIter + ".txt");
        } catch (Exception e) {
            IJ.error(e.toString());
        }
        try {
            outputStream = new PrintWriter(new FileOutputStream(results));
        } catch (FileNotFoundException e) {
            IJ.error("Could not write to results file.");
        }
//        outputStream.println("Growth Unit:," + hgu + ",Filter Radius:," + radius
//                + ",Gradient Sensitivity:," + gradSens + ",Noise:," + noise + "\n");
        outputStream.println("Iterations,Number of Tips,Total Length,Perimeter Length,Circularity,Area,Dbm,Dbs,Ds,Dss,Lacunarity,Ds R^2,Dss R^2,Dst,Dst R^2,HGU Estimate");
        Dimension dim = Toolkit.getDefaultToolkit().getScreenSize();
        ProgressDialog dialog = new ProgressDialog(null, "HGU: " + hgu + " " + thisIter + " - Growing...", false, false);
        dialog.setLocation(dim.width / 2 - dialog.getWidth() / 2, dim.height / 2 - dialog.getHeight() / 2);
        dialog.setVisible(true);
        nutrientField = new FloatProcessor(w, h);
        nutrientField.setValue(255.0);
        nutrientField.fill();
        int hyphalCount = hyphae.size();
//        ImagePlus imp = new ImagePlus("", ip);
        //imp.show();
        double hguestimate[] = new double[1000];
        for (i = 0; totalLength < maxLength; i++) {
            ip.setRoi((Roi) null);
            updateNutrientField(ip, i);
            IJ.freeMemory();
            dialog.updateProgress(totalLength, (int) maxLength);
            BatchAnalyser ba = new BatchAnalyser();
            ba.initialise(null, 0.0, 0.0, Double.MAX_VALUE, 0.0, Double.MAX_VALUE, 0.0, 1.0, 95.0);
            ba.setOutputData(BatchAnalyser.AREAS + BatchAnalyser.CIRC + BatchAnalyser.FRACTAL_DIMENSION
                    + BatchAnalyser.LACUNARITY);
            densityField = ip.toFloat(0, null);
            blurrer.blur(densityField, radius);
            h1 = h0;
            for (j = 0; j < h1; j++) {
                Hypha current = (Hypha) hyphae.get(j);
                current.grow(densityField, nutrientField, gradSens, noise);
                totalLength++;
                x = current.getX();
                y = current.getY();
                if (densityField.getPixelValue(x, y) < growThresh || nutrientField.getPixelValue(x, y) <= 0.0) {
                    hyphae.remove(j);
                    j--;
                    h0--;
                    h1--;
                } else {
                    if (x <= xlow) {
                        xlow = x - 1;
                    } else if (x >= xhigh) {
                        xhigh = x + 1;
                    }
                    if (y <= ylow) {
                        ylow = y - 1;
                    } else if (y >= yhigh) {
                        yhigh = y + 1;
                    }
                    double l = current.getLength();
//                    if (l > hgu) {
                    if (l > hgu + current.getBranchOffset()) {
                        double e = R.nextGaussian() * hgu * 0.04;
                        if (R.nextBoolean()) {
                            e *= -1;
                        }
                        apex = l * 0.2 + e;
//                        apex = l * 0.2;
                        try {
                            x = current.getBranchx(apex);
                            y = current.getBranchy(apex);
                        } catch (IndexOutOfBoundsException g) {
                            x = current.getX();
                            y = current.getY();
                        }
                        angle = current.getBranchAngle();
                        hyphae.add(new Hypha(x, y, angle, ip, hgu));
                        current.resetLength();
                        current.resetBranchOffset();
                        h0++;
                        hyphalCount++;
                    }
                }
            }
            ip.setRoi(new Rectangle(xlow, ylow, xhigh - xlow + 1, yhigh - ylow + 1));
            double dims[] = boxCounter.do2DEstimate(ip.crop());
            ba.searchImage(ip.duplicate().crop(), false, null, false);
            double params[] = ba.getParams();
            params[2] = (5.0 - Math.abs(params[2])) / 2.0;
//            hguestimate[i] = hguEstimate(params[1], params[3], params[2]);
            hguestimate[i] = 0.0;
            if (dims != null && params != null) {
//                outputStream.println(i + "," + hyphalCount + "," + totalLength + "," + decFormat.format(params[4]) + "," + decFormat.format(params[0])
//                        + "," + decFormat.format(params[1]) + "," + decFormat.format(dims[0]) + "," + decFormat.format(dims[1])
//                        + "," + decFormat.format(params[2]) + "," + decFormat.format(params[5])
//                        + "," + decFormat.format(params[3]) + "," + decFormat.format(params[6]) + "," + decFormat.format(params[7])
//                        + "," + decFormat.format(params[8]) + "," + decFormat.format(params[9]) + "," + decFormat.format(hguestimate));
                outputStream.println(i + "\t" + hyphalCount + "\t" + totalLength + "\t"
                        + decFormat.format(params[0]) + "\t" + decFormat.format(params[1])
                        + "\t" + decFormat.format(params[2]) + "\t" + decFormat.format(params[3])
                        + "\t" + decFormat.format(params[4]) + "\t" + decFormat.format(dims[0])
                        + "\t" + decFormat.format(dims[1]) + "\t" + decFormat.format(hguestimate[i]));
//            }
            } else {
//                output = i + "," + hyphalCount + ",N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A";
                outputStream.println(i + "," + hyphalCount + ",N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A");
            }
//            outputStream.println(output);
            ImageProcessor outIP = ip.duplicate();
            outIP.setFont(new Font("Arial", Font.BOLD, 30));
//            outIP.drawString(output, 10, 32);
            IJ.saveAs(new ImagePlus("", outIP), "PNG", imageFolder + "\\Mycelium" + numFormat.format(i));
//            IJ.saveAs(new ImagePlus("", densityField), "TIF", imageFolder + "\\DensityField" + numFormat.format(i));
//            IJ.saveAs(new ImagePlus("", nutrientField), "TIF", imageFolder + "\\NutrientField" + numFormat.format(i));

            if (dataPoints.size() > i) {
                double data[][] = (double[][]) dataPoints.get(i);
                int index = (int) Math.round(data[15][0]);
                data[0][index] = (double) hyphalCount;
                data[1][index] = (double) totalLength;
                if (dims != null && params != null) {
                    data[2][index] = params[4];
                    data[3][index] = params[0];
                    data[4][index] = params[1];
                    data[5][index] = dims[0];
                    data[6][index] = dims[1];
                    data[7][index] = params[2];
//                    data[8][index] = params[5];
                    data[9][index] = params[3];
//                    data[10][index] = params[6];
//                    data[11][index] = params[7];
//                    data[12][index] = params[8];
//                    data[13][index] = params[9];
                    data[14][index] = hguestimate[i];
                }
                data[15][0] += 1.0;
            } else {
                if (dims != null && params != null) {
                    double newdata[][] = new double[16][popSize];
                    newdata[0][0] = (double) hyphalCount;
                    newdata[1][0] = (double) totalLength;
                    newdata[2][0] = params[4];
                    newdata[3][0] = params[0];
                    newdata[4][0] = params[1];
                    newdata[5][0] = dims[0];
                    newdata[6][0] = dims[1];
                    newdata[7][0] = params[2];
//                    newdata[8][0] = params[5];
                    newdata[9][0] = params[3];
//                    newdata[10][0] = params[6];
//                    newdata[11][0] = params[7];
//                    newdata[12][0] = params[8];
//                    newdata[13][0] = params[9];
                    newdata[14][0] = hguestimate[i];
                    newdata[15][0] = 1.0;
                    dataPoints.add(newdata);
                } else {
                    double newdata[][] = new double[16][popSize];
                    newdata[0][0] = (double) hyphalCount;
                    newdata[1][0] = (double) totalLength;
                    newdata[2][0] = 0.0;
                    newdata[3][0] = 0.0;
                    newdata[4][0] = 0.0;
                    newdata[5][0] = 0.0;
                    newdata[6][0] = 0.0;
                    newdata[7][0] = 0.0;
                    newdata[8][0] = 0.0;
                    newdata[9][0] = 0.0;
                    newdata[10][0] = 0.0;
                    newdata[11][0] = 0.0;
                    newdata[12][0] = 0.0;
                    newdata[13][0] = 0.0;
                    newdata[14][0] = 0.0;
                    newdata[15][0] = 1.0;
                    dataPoints.add(newdata);
                }
            }
//            imp.updateAndDraw();
        }
//        System.out.print(hgu + "\t");
//        printMeanAndSD(hguestimate);
        outputStream.close();
        IJ.saveAs(new ImagePlus("", ip), "PNG", imageFolder + "\\Result_" + numFormat.format(thisIter));
//        imp.close();
        dialog.dispose();
    }

//    void printMeanAndSD(double data[]) {
//        double sum = 0.0d;
//        double mean, sd;
//        int count = 0;
//        for (int i = 0; i < data.length; i++) {
//            if (!(Double.isNaN(data[i])) && data[i] > 0.0) {
//                sum += data[i];
//                count++;
//            }
//        }
//        mean = sum / count;
//        double sumvar = 0.0d;
//        for (int j = 0; j < data.length; j++) {
//            if (!(Double.isNaN(data[j])) && data[j] > 0.0) {
//                sumvar += Math.pow(data[j] - mean, 2);
//            }
//        }
//        double variance = sumvar / count;
//        sd = Math.sqrt(variance);
//        System.out.print(mean + "\t" + sd + "\n");
//    }
    void updateNutrientField(ByteProcessor ip, int i) {
        FloatBlitter nutFieldBlit = new FloatBlitter(nutrientField);
        nutFieldBlit.copyBits(ip, 0, 0, FloatBlitter.MIN);
        (new GaussianBlur()).blur(nutrientField, radius * 4.0);
    }

    static double[] getMeanAndSD(double data[]) {
        double result[] = new double[2];
        double sum = 0.0;
        for (int i = 0; i < data.length; i++) {
            sum += data[i];
        }
        result[0] = sum / data.length;
        sum = 0.0;
        for (int i = 0; i < data.length; i++) {
            sum += Math.pow(data[i] - result[0], 2.0);
        }
        result[1] = Math.sqrt(sum / data.length);
        return result;
    }

    double hguEstimate(double area, double lac, double ds) {
        double areaest = Math.log(area) / area_max;
        int areaindex = 0;
        while (areaindex < areaCurve.length && areaCurve[areaindex] < areaest) {
            areaindex++;
        }
        if (areaindex < 10) {
            areaindex = 10;
        } else if (areaindex >= areaCurve.length - 10) {
            areaindex = areaCurve.length - 11;
        }
        double vector1[] = new double[2];
        double vector2[] = new double[2];
        vector1[0] = lac / lac_max;
        vector1[1] = ds / ds_max;
//        vector1[2] = Math.log(area) / area_max;
        double minDist = Double.MAX_VALUE;
        int hguindex = -1;
        for (int i = 0; i < lacSurf[0].getWidth(); i++) {
            for (int j = areaindex - 10; j <= areaindex + 10; j++) {
//            for (int j = 0; j < lacSurf[0].getHeight(); j++) {
                int imageindex = R.nextInt(10);
                vector2[0] = lacSurf[imageindex].getPixelValue(i, j);
                imageindex = R.nextInt(10);
                vector2[1] = dsSurf[imageindex].getPixelValue(i, j);
//                imageindex = R.nextInt(10);
//                vector2[2] = areaSurf[imageindex].getPixelValue(i, j);
                double dist = Utils.calcEuclidDist(vector1, vector2);
                if (dist < minDist) {
                    minDist = dist;
                    hguindex = i;
                }
            }
        }
        if (hguindex > -1) {
            return 10.0 + 5.0 * hguindex;
        } else {
            return Double.NaN;
        }
    }
}
