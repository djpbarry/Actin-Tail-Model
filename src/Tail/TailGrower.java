package Tail;

import UtilClasses.Utilities;
import IAClasses.ProgressDialog;
import IAClasses.Utils;
import UtilClasses.GenUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.process.ByteProcessor;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;

public class TailGrower {

    public double radius = 25.0;
//    private FloatProcessor densityField, nutrientField;
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
    private double res = 7; // One pixel equals 7 nm, approximate width of actin filament - http://www.ncbi.nlm.nih.gov/books/NBK9908/
    private double virusRadius = 250.0 / res;
    private double capFactor = 1.0;

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
        Random rand = new Random();
        double apex;
        DecimalFormat numFormat = new DecimalFormat("000");
        DecimalFormat decFormat = new DecimalFormat("0.000");
        int totalLength = 0;
        ip.setValue(0);
        ArrayList<Filament> filaments = new ArrayList();
        Virus virus = new Virus(ip.getWidth() / 2.0, ip.getHeight() / 2.0);
        filaments.add(createInitialFilament(virus, ip, branchRate));
        int filCount = 1;
        imageFolder = GenUtils.createDirectory(((Utilities.getFolder(imageFolder, "Choose_Location_for_Output", true)).getAbsolutePath()
                + "\\Sim_Actin_Tails\\" + decFormat.format(branchRate) + "_" + maxLength));
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
        int fCount = filaments.size();
        double logMax = Math.log(maxLength);
        int virRadius = (int) Math.round(virus.getRadius() / res);
        for (int i = 0; totalLength < maxLength && filaments.size() > 0; i++) {
//            ip.setRoi((Roi) null);
            ip.drawOval((int) Math.round(virus.getX() / res) - virRadius,
                    (int) Math.round(virus.getY() / res) - virRadius,
                    2 * virRadius, 2 * virRadius);
            IJ.freeMemory();
            dialog.updateProgress(Math.log(totalLength), logMax);
            int f1 = filCount;
            Force totalForce = new Force(0.0, 0.0);
            for (int j = 0; j < f1; j++) {
                Filament current = filaments.get(j);
                int x = (int) Math.round(current.getX() / res);
                int y = (int) Math.round(current.getY() / res);
                if (current.grow(noise, virus)) {
                    totalLength++;
                    x = (int) Math.round(current.getX() / res);
                    y = (int) Math.round(current.getY() / res);
                    ip.drawPixel(x, y);
                }
                double clearance = Utils.calcDistance(x, y, virus.getX(), virus.getY());
                if (clearance > rand.nextDouble() * res * branchRate * capFactor) {
                    filaments.remove(j);
                    j--;
                    filCount--;
                    f1--;
                    if (filCount < 1) {
                        filaments.add(createInitialFilament(virus, ip, branchRate));
                        filCount++;
                    }
                } else {
                    totalForce.addForce(current.calcForce(virus.getX(), virus.getY(), virus.getRadius()));
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
                        double angle = current.getBranchAngle();
                        filaments.add(new Filament(x, y, angle, branchRate, res));
                        current.resetLength();
                        current.resetBranchOffset();
                        filCount++;
                        fCount++;
                    }
                }
            }
            virus.updateVelocity(totalForce);
            virus.updatePosition();
//            ip.setRoi(new Rectangle(xlow, ylow, xhigh - xlow + 1, yhigh - ylow + 1));
            outputStream.println(i + "\t" + fCount + "\t" + totalLength);
            ImageProcessor outIP = ip.duplicate();
            if (showAllImages) {
                IJ.saveAs(new ImagePlus("", outIP), "PNG", imageFolder + "\\Step" + numFormat.format(i));
            }
        }
        outputStream.close();
        IJ.saveAs(new ImagePlus("", ip), "PNG", imageFolder + "\\Result");
        dialog.dispose();
    }

    Filament createInitialFilament(Virus virus, ImageProcessor ip, double branchRate) {
        Random rand = new Random();
        double angle = 2.0 * Math.PI * rand.nextDouble();
        double r = Math.abs(rand.nextGaussian()) * virus.getRadius() * 0.01 + virus.getRadius();
        double x0 = (ip.getWidth() / 2.0 + r * Math.cos(angle)) * res;
        double y0 = (ip.getHeight() / 2.0 + r * Math.sin(angle)) * res;
        return new Filament(x0, y0, 360.0 * rand.nextDouble(), branchRate, res);
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
