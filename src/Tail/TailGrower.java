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
import org.apache.commons.math3.analysis.function.Gaussian;

public class TailGrower {

//    private FloatProcessor densityField, nutrientField;
    public double noise = 20.0;
    public File imageFolder = new File("c:/users/barry05/desktop");
    public double areaCurve[];
    public double ds_max = -Double.MAX_VALUE, lac_max = -Double.MAX_VALUE,
            area_max = -Double.MAX_VALUE;
    private Random R = new Random();
    private static int maxSteps = 1500,
            width = 1500, height = 1500;
    private final String TITLE = "",
            resultsHeadings = "Iterations\tNumber of Tips\tTotal Length";
    private static boolean showAllImages = true;
    private double res = 7; // One pixel equals 7 nm, approximate width of actin filament - http://www.ncbi.nlm.nih.gov/books/NBK9908/
    private double capFactor = 150.0;
    private int simultaneousFils = 20;
    private double pZoneDegrees = 5.0;
    private double minFilLengthForBranch = 50.0;
    private double branchZoneWidth = res;
    private double t = 1.0;
    Random rand = new Random();

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
        grow(bp, maxSteps);
    }

    public TailGrower() {
    }

    public void grow(ByteProcessor ip, int maxSteps) {
        DecimalFormat numFormat = new DecimalFormat("000");
        int totalLength = 0;
        ip.setValue(0);
        ArrayList<Filament> filaments = new ArrayList();
        Virus virus = new Virus(res * ip.getWidth() / 2.0, res * ip.getHeight() / 2.0);
        for (int i = 0; i < simultaneousFils; i++) {
            filaments.add(createInitialFilament(virus, ip));
        }
        int filCount = 1;
        imageFolder = new File(GenUtils.openResultsDirectory((Utilities.getFolder(imageFolder, "Choose_Location_for_Output", true)).getAbsolutePath()
                + "\\Sim_Actin_Tails\\" + "_" + maxSteps, "\\"));
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
        ProgressDialog dialog = new ProgressDialog(null, "Growing...", false, TITLE, false);
        dialog.setVisible(true);
        int fCount = filaments.size();
        int virRadius = (int) Math.round(virus.getRadius() / res);
        for (int i = 0; i < maxSteps; i++) {
            virus.brownian();
            IJ.freeMemory();
            dialog.updateProgress(i, maxSteps);
            int f1 = filCount;
            Energy netEnergy = new Energy(0.0, 0.0);
            for (int j = 0; j < f1; j++) {
                Filament current = filaments.get(j);
                if (current.grow(noise, virus)) {
                    totalLength++;
                    int x = (int) Math.round(current.getX() / res);
                    int y = (int) Math.round(current.getY() / res);
                    ip.drawPixel(x, y);
                }
                double dist = Utils.calcDistance(current.getX(), current.getY(),
                        virus.getX(), virus.getY()) - virus.getRadius();
                Energy currentPE = current.calcPE(virus.getX(), virus.getY(), virus.getRadius());
                if (dist > ((rand.nextDouble() * branchZoneWidth * capFactor / res)
                        * getGrowP(virus, current.getX(), current.getY()))
                        || currentPE.getMag() > Filament.MAX_FIL_PE) {
                    filaments.remove(j);
                    j--;
                    filCount--;
                    f1--;
                    if (filCount < simultaneousFils) {
                        filaments.add(createInitialFilament(virus, ip));
                        filCount++;
                    }
                } else {
                    netEnergy.addEnergy(current.calcPE(virus.getX(), virus.getY(), virus.getRadius()));
                    double l = current.getLength() * current.getThickness();
                    if (l > minFilLengthForBranch && rand.nextDouble() > getBranchP(dist)) {
                        double x = current.getX();
                        double y = current.getY();
                        double angle = current.getBranchAngle();
                        filaments.add(new Filament(x, y, angle));
                        current.resetLength();
                        filCount++;
                        fCount++;
                    }
                }
            }
            virus.updateVelocity(netEnergy, t);
            outputStream.println(i + "\t" + fCount + "\t" + totalLength);
            ImageProcessor filOut = ip.duplicate();
            ImageProcessor virOut = ip.duplicate();
            virOut.setValue(255);
            virOut.fill();
            virOut.setValue(0);
            virOut.fillOval((int) Math.round(virus.getX() / res) - virRadius,
                    (int) Math.round(virus.getY() / res) - virRadius,
                    2 * virRadius, 2 * virRadius);
//            ByteProcessor growthZone = new ByteProcessor(ip.getWidth(), ip.getHeight());
//            growthZone.setValue(0);
//            growthZone.fill();
//            double x0 = virus.getX();
//            double y0 = virus.getY();
//            double r = 2.0 * virus.getRadius();
//            for (double j = y0 - r; j <= y0 + r; j += res) {
//                for (double k = x0 - r; k <= x0 + r; k += res) {
//                    growthZone.putPixelValue((int) Math.round(k / res),
//                            (int) Math.round(j / res), (int) Math.round(255.0 * getGrowP(virus, k, j)));
//                }
//            }
            if (showAllImages) {
                IJ.saveAs(new ImagePlus("", filOut), "PNG", imageFolder + "\\Filaments_Step" + numFormat.format(i));
                IJ.saveAs(new ImagePlus("", virOut), "PNG", imageFolder + "\\Virus_Step" + numFormat.format(i));
//                IJ.saveAs(new ImagePlus("", growthZone), "PNG", imageFolder + "\\GrowthZone_Step" + numFormat.format(i));
            }
        }
        outputStream.close();
        IJ.saveAs(new ImagePlus("", ip), "PNG", imageFolder + "\\Result");
        dialog.dispose();
    }

    Filament createInitialFilament(Virus virus, ImageProcessor ip) {
        double angle = Math.toRadians(getInitAngle(virus));
        double r = rand.nextDouble() * branchZoneWidth / res + virus.getRadius();
        double x0 = virus.getX() + r * Math.cos(angle);
        double y0 = virus.getY() - r * Math.sin(angle);
        return new Filament(x0, y0, 360.0 * rand.nextDouble());
    }

    private double getInitAngle(Virus virus) {
        double avVel[] = virus.getVelAv(100);
        double angle = Utils.arcTan(avVel[0], avVel[1]) - 180.0 + pZoneDegrees * rand.nextGaussian();
        if (angle < 0.0) {
            angle += 360.0;
        }
        return angle;
    }

    private double getGrowP(Virus virus, double x, double y) {
        double avVel[] = virus.getVelAv(100);
        double peakPos = Utils.arcTan(avVel[0], avVel[1]) + 180.0;
        if (peakPos > 360.0) {
            peakPos -= 360.0;
        }
        double pos = Utils.arcTan(x - virus.getX(), y - virus.getY());
        if (pos > 360.0) {
            pos -= 360.0;
        }
        if (pos >= 0.0 && pos < 90.0 && peakPos >= 270.0) {
            peakPos -= 180.0;
            pos += 180.0;
        } else if (peakPos >= 0.0 && peakPos < 90.0 && pos >= 270.0) {
            pos -= 180.0;
            peakPos += 180.0;
        }
        Gaussian gaussian = new Gaussian(0.0, pZoneDegrees);
        return gaussian.value(pos - peakPos) / gaussian.value(0.0);
    }

    private double getBranchP(double d) {
        Gaussian gaussian = new Gaussian(0.0, branchZoneWidth);
        return gaussian.value(d) / gaussian.value(0.0);
    }

    public boolean showDialog(String title) {
        boolean valid = false;
        while (!valid) {
            GenericDialog dialog = new GenericDialog(title, IJ.getInstance());
            dialog.addNumericField("Maximum length:", maxSteps, 0, 5, "Pixels");
            dialog.addNumericField("Image width:", width, 0, 5, "");
            dialog.addNumericField("Image height:", height, 0, 5, "");
            dialog.addCheckbox("Show all images?", showAllImages);
            if (IJ.getInstance() != null) {
                dialog.setLocation((IJ.getInstance()).getPreferredLocation());
            }
            dialog.showDialog();
            if (!dialog.wasCanceled()) {
                maxSteps = (int) Math.round(dialog.getNextNumber());
                width = (int) Math.round(dialog.getNextNumber());
                height = (int) Math.round(dialog.getNextNumber());
                showAllImages = dialog.getNextBoolean();
                if (dialog.invalidNumber()) {
                    GenUtils.error("Entries must be numeric.");
                } else if (width <= 0 || height <= 0) {
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
