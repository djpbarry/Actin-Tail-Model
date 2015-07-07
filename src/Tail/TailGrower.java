package Tail;

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
import ij.process.TypeConverter;
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

public class TailGrower {

    public double radius = 25.0;
    private FloatProcessor densityField, nutrientField;
    public double dfGradSens = 5000.0, nfGradSens = 500.0, noise = 1.5, growThresh = 150.0;
    public File imageFolder = new File("c:/users/barry05/desktop");
    public ImageProcessor[] lacSurf, dsSurf, areaSurf;
    public double areaCurve[];
    public double ds_max = -Double.MAX_VALUE, lac_max = -Double.MAX_VALUE,
            area_max = -Double.MAX_VALUE;
    private Random R = new Random();
    private static int hgu = 20, maxLength = 400000,
            width = 1000, height = 500;
    private final String TITLE = "",
            resultsHeadings = "Iterations\tNumber of Tips\tTotal Length";
    private static boolean showAllImages = true;
    private double res = 7; // One pixel equals 7 nm, approximate width of actin filament
    private double virusRadius = 250.0/res;

    public static void main(String args[]) {
        TailGrower grower = new TailGrower();
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
        ByteProcessor bp = new ByteProcessor(width, height);
        bp.setColor(Color.white);
        bp.fill();
        grow(bp, hgu, maxLength);
    }

    public TailGrower() {
    }

    public void grow(ByteProcessor ip, double branchRate, int maxLength) {
        int w = ip.getWidth(), h = ip.getHeight(), i, j, h0 = 0, h1, x0 = 2,
                y0 = h / 2, x, y, xlow = w / 2 - 1, xhigh = w / 2 + 1,
                ylow = h / 2 - 1, yhigh = h / 2 + 1;
        double apex;
        DecimalFormat numFormat = new DecimalFormat("000");
        DecimalFormat decFormat = new DecimalFormat("0.000");
        int totalLength = 0;

        ip.setValue(0);
        ArrayList filaments = new ArrayList();
        double angle = 0.0;

        filaments.add(new Filament(x0, y0, angle, ip, branchRate));
        h0++;
        imageFolder = GenUtils.createDirectory(((Utilities.getFolder(imageFolder, "Choose_Location_for_Output", true)).getAbsolutePath()
                + "\\Mycelium\\" + decFormat.format(branchRate) + "_" + maxLength));
        File results = null;
        PrintWriter outputStream = null;
        try {
            results = new File(imageFolder + "\\Results.txt");
        } catch (Exception e) {
            IJ.error(e.toString());
        }
        try {
            outputStream = new PrintWriter(new FileOutputStream(results));
        } catch (FileNotFoundException e) {
            IJ.error("Could not write to results file.");
        }
        outputStream.println(resultsHeadings);
        ProgressDialog dialog = new ProgressDialog(null, "Branch Rate: " + branchRate + " - Growing...", false, TITLE, false);
        dialog.setVisible(true);
        nutrientField = new FloatProcessor(w, h);
        for (x = 0; x < w; x++) {
            for (y = 0; y < h; y++) {
                nutrientField.putPixelValue(x, y, 1.0 - Math.abs(0.5 * h - y) / (0.5 * h));
            }
        }
        IJ.saveAs(new ImagePlus("", nutrientField), "TIF", imageFolder + "\\NF_Step_Init");
        int fCount = filaments.size();
        double logMax = Math.log(maxLength);
        double xmax = x0;
        for (i = 0; totalLength < maxLength && filaments.size() > 0; i++) {
            ip.setRoi((Roi) null);
//            updateNutrientField(ip, i);
            IJ.freeMemory();
            dialog.updateProgress(Math.log(totalLength), logMax);
//            densityField = ip.toFloat(0, null);
//            blurrer.blur(densityField, radius);
            h1 = h0;
            for (j = 0; j < h1; j++) {
                Filament current = (Filament) filaments.get(j);
                current.grow(densityField, nutrientField, dfGradSens, noise, nfGradSens);
                totalLength++;
                x = current.getX();
                y = current.getY();
                if (x > xmax) {
                    xmax = x;
                }
//                if (densityField.getPixelValue(x, y) < growThresh || nutrientField.getPixelValue(x, y) <= 0.0) {
                if (Math.abs(y - 0.5 * h) / h > 0.2 + 0.04 * R.nextGaussian()
                        || (xmax - x) > (0.1 + 0.02 * R.nextGaussian()) * w) {
                    filaments.remove(j);
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
                    if (l > branchRate + current.getBranchOffset()) {
                        double e = R.nextGaussian() * branchRate * 0.04;
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
                        filaments.add(new Filament(x, y, angle, ip, branchRate));
                        current.resetLength();
                        current.resetBranchOffset();
                        h0++;
                        fCount++;
                    }
                }
            }
            ip.setRoi(new Rectangle(xlow, ylow, xhigh - xlow + 1, yhigh - ylow + 1));
            outputStream.println(i + "\t" + fCount + "\t" + totalLength);
            ImageProcessor outIP = ip.duplicate();
            outIP.setFont(new Font("Arial", Font.BOLD, 30));
            if (showAllImages) {
                IJ.saveAs(new ImagePlus("", outIP), "PNG", imageFolder + "\\Step" + numFormat.format(i));
//                IJ.saveAs(new ImagePlus("", nutrientField), "TIF", imageFolder + "\\NF_Step" + numFormat.format(i));
            }
        }
        outputStream.close();
        IJ.saveAs(new ImagePlus("", ip), "PNG", imageFolder + "\\Result");
        dialog.dispose();
    }

    void updateNutrientField(ByteProcessor ip, int i) {
        FloatBlitter nutFieldBlit = new FloatBlitter(nutrientField);
        FloatProcessor ipDup = (FloatProcessor) (new TypeConverter(ip, false)).convertToFloat(null);
        ipDup.invert();
        ipDup.multiply(1.0 / 255.0);
        nutFieldBlit.copyBits(ipDup, 0, 0, FloatBlitter.SUBTRACT);
        (new GaussianBlur()).blur(nutrientField, radius);
    }

    public boolean showDialog(String title) {
        boolean valid = false;
        while (!valid) {
            GenericDialog dialog = new GenericDialog(title, IJ.getInstance());
            dialog.addNumericField("Branching Constant:", hgu, 0, 5, "");
            dialog.addNumericField("Maximum length:", maxLength, 0, 5, "Pixels");
            dialog.addNumericField("Image width:", width, 0, 5, "");
            dialog.addNumericField("Image height:", height, 0, 5, "");
            dialog.addCheckbox("Show all images?", showAllImages);
            if (IJ.getInstance() != null) {
                dialog.setLocation((IJ.getInstance()).getPreferredLocation());
            }
            dialog.showDialog();
            if (!dialog.wasCanceled()) {
                hgu = (int) Math.round(dialog.getNextNumber());
                maxLength = (int) Math.round(dialog.getNextNumber());
                width = (int) Math.round(dialog.getNextNumber());
                height = (int) Math.round(dialog.getNextNumber());
                showAllImages = dialog.getNextBoolean();
                if (dialog.invalidNumber()) {
                    GenUtils.error("Entries must be numeric.");
                } else if (width <= 0 || height <= 0 || hgu < 0) {
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
