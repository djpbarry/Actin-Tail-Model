package Mycelium;

import UtilClasses.Utilities;
import IAClasses.ProgressDialog;
import UtilClasses.GenUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.Roi;
import ij.plugin.filter.GaussianBlur;
import ij.process.ByteProcessor;
import ij.process.FloatBlitter;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.awt.Rectangle;
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
    private static int firstH = 50, lastH = 50, intervalH = 5, maxLength = 40000,
            width = 1200, height = 1200;
    private final String TITLE = "",
            resultsHeadings = "Iterations\tNumber of Tips\tTotal Length";
    private static boolean showAllImages = false;

    public static void main(String args[]) {
        MyceliumGrower grower = new MyceliumGrower();
        if (!grower.initialise()) {
            System.exit(-1);
        }
        grower.run();
        System.exit(0);
    }

    boolean initialise() {
        if (!showDialog(TITLE)) {
            return false;
        }
        return true;
    }

    void run() {
        for (int h = firstH; h <= lastH; h += intervalH) {
            dataPoints = new ArrayList();
            for (int i = 0; i < popSize; i++) {
                ByteProcessor bp = new ByteProcessor(width, height);
                bp.setColor(Color.white);
                bp.fill();
                growMycelium(bp, i, h, maxLength);
            }
        }
    }

    public MyceliumGrower() {
    }

    public void growMycelium(ByteProcessor ip, int thisIter, double hgu, int maxLength) {
        int w = ip.getWidth(), h = ip.getHeight(), i, j, h0 = 0, h1, x0 = w / 2,
                y0 = h / 2, x, y, xlow = w / 2 - 1, xhigh = w / 2 + 1,
                ylow = h / 2 - 1, yhigh = h / 2 + 1;
        double apex;
        GaussianBlur blurrer = new GaussianBlur();
        DecimalFormat numFormat = new DecimalFormat("000");
        DecimalFormat decFormat = new DecimalFormat("0.000");
        int totalLength = 0;

        ip.setValue(0);
        ArrayList hyphae = new ArrayList();
        double angle = R.nextDouble() * 360;

        hyphae.add(new Hypha(x0, y0, angle, ip, hgu));
        h0++;
        imageFolder = GenUtils.createDirectory(((Utilities.getFolder(imageFolder, "Choose_Location_for_Output")).getAbsolutePath()
                + "\\Mycelium\\" + decFormat.format(hgu) + "_" + maxLength));
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
        outputStream.println(resultsHeadings);
        ProgressDialog dialog = new ProgressDialog(null, "HGU: " + hgu + " " + thisIter + " - Growing...", false, false, TITLE);
        dialog.setVisible(true);
        nutrientField = new FloatProcessor(w, h);
        nutrientField.setValue(255.0);
        nutrientField.fill();
        int hyphalCount = hyphae.size();
        double logMax = Math.log(maxLength);
        for (i = 0; totalLength < maxLength; i++) {
            ip.setRoi((Roi) null);
            updateNutrientField(ip, i);
            IJ.freeMemory();
            dialog.updateProgress(Math.log(totalLength), logMax);
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
            outputStream.println(i + "\t" + hyphalCount + "\t" + totalLength);
            ImageProcessor outIP = ip.duplicate();
            outIP.setFont(new Font("Arial", Font.BOLD, 30));
            if (showAllImages) {
                IJ.saveAs(new ImagePlus("", outIP), "PNG", imageFolder + "\\Step" + numFormat.format(i));
            }
        }
        outputStream.close();
        IJ.saveAs(new ImagePlus("", ip), "PNG", imageFolder + "\\Result_" + numFormat.format(thisIter));
        dialog.dispose();
    }

    void updateNutrientField(ByteProcessor ip, int i) {
        FloatBlitter nutFieldBlit = new FloatBlitter(nutrientField);
        nutFieldBlit.copyBits(ip, 0, 0, FloatBlitter.MIN);
        (new GaussianBlur()).blur(nutrientField, radius * 4.0);
    }

    public boolean showDialog(String title) {
        boolean valid = false;
        while (!valid) {
            GenericDialog dialog = new GenericDialog(title, IJ.getInstance());
            dialog.addNumericField("Begin at:", firstH, 0, 5, "");
            dialog.addNumericField("End at:", lastH, 0, 5, "");
            dialog.addNumericField("Step size:", intervalH, 0, 5, "");
            dialog.addNumericField("Maximum length:", maxLength, 0, 5, "Pixels");
            dialog.addNumericField("Image width:", width, 0, 5, "");
            dialog.addNumericField("Image height:", height, 0, 5, "");
            dialog.addCheckbox("Show all images?", showAllImages);
            if (IJ.getInstance() != null) {
                dialog.setLocation((IJ.getInstance()).getPreferredLocation());
            }
            dialog.showDialog();
            if (!dialog.wasCanceled()) {
                firstH = (int) Math.round(dialog.getNextNumber());
                lastH = (int) Math.round(dialog.getNextNumber());
                intervalH = (int) Math.round(dialog.getNextNumber());
                maxLength = (int) Math.round(dialog.getNextNumber());
                width = (int) Math.round(dialog.getNextNumber());
                height = (int) Math.round(dialog.getNextNumber());
                showAllImages = dialog.getNextBoolean();
                if (dialog.invalidNumber()) {
                    GenUtils.error("Entries must be numeric.");
                } else if (width <= 0 || height <= 0 || firstH < 0 || lastH < 0 || intervalH <= 0) {
                    GenUtils.error("Values must be greater than zero.");
                } else {
                    valid = true;
                }
            } else {
                return false;
            }
        }
        return true;
    }
}
