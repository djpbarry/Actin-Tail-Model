package Mycelium;

import AnaMorf.BatchAnalyser;
import IAClasses.FractalEstimator;
import EMSeg.ProgressDialog;
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
    private static int popSize = 10;
    public File imageFolder;
    public ImageProcessor lacSurf, fracSurf, areaCurve;
    private final double DS_MAX = 1.38, LAC_MAX = 6.543;

    public static void main(String args[]) {
        MyceliumGrower grower = new MyceliumGrower();
        float steps = 40000.0f;
        Random rand = new Random();
        double maxinc = 10.0;
        grower.lacSurf = (new ImagePlus("C:\\Users\\Dave\\Desktop\\lac.tif")).getProcessor();
        grower.fracSurf = (new ImagePlus("C:\\Users\\Dave\\Desktop\\ds.tif")).getProcessor();
        grower.areaCurve = (new ImagePlus("C:\\Users\\Dave\\Desktop\\area.tif")).getProcessor();
        for (int x = 0; x < grower.areaCurve.getWidth(); x++) {
            double sum = 0.0;
            for (int y = 0; y < grower.areaCurve.getHeight(); y++) {
                sum += grower.areaCurve.getPixelValue(x, y);
            }
            grower.areaCurve.putPixelValue(x, 1, sum / grower.areaCurve.getHeight());
        }
        for (int h = 85; h <= 100; h += 5) {
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
        //ImagePlus imp = new ImagePlus("", ip);
        int w = ip.getWidth(), h = ip.getHeight(), i, j, h0 = 0, h1, x0 = w / 2,
                y0 = h / 2, x, y, xlow = w / 2 - 1, xhigh = w / 2 + 1,
                ylow = h / 2 - 1, yhigh = h / 2 + 1;
        double apex;
        Random R = new Random();
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

        imageFolder = new File("d:\\Mycelium\\" + decFormat.format(hgu)
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
//        outputStream.println("Iterations,Number of Tips,Total Length,Perimeter Length,Circularity,Area,Dbm,Dbs,Ds,Dss,Lacunarity,Ds R^2,Dss R^2,Dst,Dst R^2,HGU Estimate");
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
                current.grow();
                totalLength++;
                x = current.getX();
                y = current.getY();
                /*
                 * if (densityField.getPixelValue(x, y) < growThresh ||
                 * nutrientField.getPixelValue(x, y) <= 0.0) { hyphae.remove(j);
                 * j--; h0--; h1--; } else
                 */ {
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
                    if (l > hgu + current.getBranchOffset()) {
                        double e = R.nextGaussian() * hgu * 0.04;
                        if (R.nextBoolean()) {
                            e *= -1;
                        }
                        apex = l * 0.2 + e;
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
//            double hguestimate = hguEstimate(params[1], params[3], params[2]);
            double hguestimate = 0.0;
//            System.out.println(hguestimate);
//            double hguestimate = 0.0;
//            String output;
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
                        + "\t" + decFormat.format(dims[1]) + "\t" + decFormat.format(hguestimate));
            }
//            } else {
//                output = i + "," + hyphalCount + ",N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A";
//            }
//            outputStream.println(output);
            //ImageProcessor outIP = ip.duplicate();
            //outIP.setFont(new Font("Arial", Font.BOLD, 30));
            //outIP.drawString(output, 10, 32);
            //IJ.saveAs(new ImagePlus("", outIP), "PNG", imageFolder + "\\Mycelium" + numFormat.format(i));
            //IJ.saveAs(new ImagePlus("", densityField), "TIF", imageFolder + "\\DensityField" + numFormat.format(i));
            //IJ.saveAs(new ImagePlus("", nutrientField), "TIF", imageFolder + "\\NutrientField" + numFormat.format(i));
            //System.out.println("x: "+totalLength+" N: "+hyphalCount+ " HGU: "+(float)totalLength/hyphalCount);

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
                    data[14][index] = hguestimate;
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
                    newdata[14][0] = hguestimate;
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
            //imp.updateAndDraw();
        }
        outputStream.close();
        IJ.saveAs(new ImagePlus("", ip), "PNG", imageFolder + "\\Result_" + numFormat.format(thisIter));
//        imp.close();
        dialog.dispose();
        return;
    }

    void updateNutrientField(ByteProcessor ip, int i) {
        FloatBlitter nutFieldBlit = new FloatBlitter(nutrientField);
        nutFieldBlit.copyBits(ip, 0, 0, FloatBlitter.MIN);
        (new GaussianBlur()).blur(nutrientField, radius * 4.0);
        return;
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

    double hguEstimate(double area, double lac, double frac) {
//        int aIndex = 0;
//        while (areaCurve.getPixelValue(aIndex, 1) < area && aIndex < areaCurve.getWidth()) {
//            aIndex++;
//        }
        double vector1[] = new double[2];
        double vector2[] = new double[2];
//        vector1[0] = (lac - LAC_MIN) / LAC_MAX;
//        vector1[1] = (frac - DS_MIN) / DS_MAX;
        vector1[0] = lac / LAC_MAX;
        vector1[1] = frac / DS_MAX;
        double minDist = Double.MAX_VALUE;
        int hguindex = -1;
        for (int i = 0; i < lacSurf.getWidth(); i++) {
            for (int j = 0; j < lacSurf.getHeight(); j++) {
                vector2[0] = lacSurf.getPixelValue(i, j);
                vector2[1] = fracSurf.getPixelValue(i, j);
                double dist = Utils.calcEuclidDist(vector1, vector2);
                if (dist < minDist) {
                    minDist = dist;
                    hguindex = j;
                }
            }
        }
        if (hguindex > -1) {
            return 100.0 - 5.0 * hguindex;
        } else {
            return Double.NaN;
        }
    }

    public class Hypha {

        private double x, y, hgu;
        int length;
        ArrayList xPix = new ArrayList();
        ArrayList yPix = new ArrayList();
        ImageProcessor plot;
        Random R = new Random();
        boolean branchX = true;
        private double angle, branch = -90, branchOffset;

        public Hypha(double xc, double yc, double a0, ImageProcessor ip, double hgu) {
            this.x = xc;
            this.y = yc;
            this.angle = a0;
            this.plot = ip;
            this.length = 0;
            this.hgu = hgu;
            branchOffset = R.nextGaussian() * hgu * 0.2;
            if (R.nextBoolean()) {
                branchOffset *= -1;
            }
        }

        public void grow() {
            xPix.add(new Double(x));
            yPix.add(new Double(y));
            double nfParams[] = getVector(nutrientField, x, y);
            double dfParams[] = getVector(densityField, x, y);
            if (nfParams == null) {
                nfParams = new double[2];
                nfParams[0] = 0.0;
                nfParams[1] = angle;
            }
            if (dfParams == null) {
                dfParams = new double[2];
                dfParams[0] = 0.0;
                dfParams[1] = angle;
            }
            angle += ((dfParams[0] * (dfParams[1] - angle) + nfParams[0]
                    * (nfParams[1] - angle)) / gradSens) + noise * R.nextGaussian();

            double xVec = Math.cos(Math.toRadians(angle));
            double yVec = Math.sin(Math.toRadians(angle));

            x += xVec;
            y += yVec;
            plot.drawLine((int) Math.round(((Double) xPix.get(length)).doubleValue()),
                    (int) Math.round(((Double) yPix.get(length)).doubleValue()),
                    (int) Math.round(x), (int) Math.round(y));
            length++;
        }

        public int getLength() {
            return length;
        }

        public void resetLength() {
            Double tempX = ((Double) xPix.get(length - 1));
            Double tempY = ((Double) yPix.get(length - 1));
            xPix.clear();
            yPix.clear();
            xPix.add(tempX);
            yPix.add(tempY);
            length = 1;
        }

        public int getBranchx(double d) {
            return (int) Math.round((Double) xPix.get(length - (int) Math.round(d) - 1));
        }

        public int getBranchy(double d) {
            return (int) Math.round((Double) yPix.get(length - (int) Math.round(d) - 1));
        }

        public double getBranchAngle() {
            branch += -2 * branch;
            return (branch + angle);
        }

        public int getX() {
            return (int) x;
        }

        public int getY() {
            return (int) y;
        }

        public double getBranchOffset() {
            return branchOffset;
        }

        public void resetBranchOffset() {
            branchOffset = R.nextGaussian() * hgu * 0.2;
            if (R.nextBoolean()) {
                branchOffset *= -1;
            }
        }

        double[] getVector(FloatProcessor field, double x, double y) {
            double xGrad = field.getInterpolatedValue(x + 1.0, y - 1.0)
                    + 2.0 * field.getInterpolatedValue(x + 1.0, y)
                    + field.getInterpolatedValue(x + 1.0, y + 1.0)
                    - field.getInterpolatedValue(x - 1.0, y - 1.0)
                    - 2.0 * field.getInterpolatedValue(x - 1.0, y)
                    - field.getInterpolatedValue(x - 1.0, y + 1.0);
            double yGrad = field.getInterpolatedValue(x - 1.0, y + 1.0)
                    + 2.0 * field.getInterpolatedValue(x, y + 1.0)
                    + field.getInterpolatedValue(x + 1.0, y + 1.0)
                    - field.getInterpolatedValue(x - 1.0, y - 1.0)
                    - 2.0 * field.getInterpolatedValue(x, y - 1.0)
                    - field.getInterpolatedValue(x + 1.0, y - 1.0);
            double gradMag = Math.sqrt(xGrad * xGrad + yGrad * yGrad);

            double t;
            if (gradMag > 0.0) {
                t = Math.toDegrees(Math.atan(yGrad / xGrad));
            } else {
                t = Double.NaN;
            }
            if (xGrad > 0.0 && yGrad < 0.0) {
                t += 360.0;
            } else if (xGrad < 0.0) {
                t += 180.0;
            }
            double params[] = {gradMag, t};
            if (Double.isNaN(t)) {
                params = null;
            }

            return params;
        }
    }
}
