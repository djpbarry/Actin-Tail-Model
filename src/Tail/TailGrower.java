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
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;
import org.apache.commons.math3.analysis.function.Gaussian;

public class TailGrower {

//    private FloatProcessor densityField, nutrientField;
    private final double FIL_NOISE = 20.0;
    public File imageFolder = new File("c:/users/barry05/desktop");
    public double areaCurve[];
    private final Random R = new Random();
    private static int maxSteps = 1000,
            width = 1000, height = 1000;
    private final String TITLE = "";
    private static boolean showAllImages = true;
    private final double RES = 7; // One pixel equals 7 nm, approximate width of actin filament - http://www.ncbi.nlm.nih.gov/books/NBK9908/
    private final double CAP_FAC = 150.0;
    private final int N_FILS = 20;
    private final double P_ZONE_DEG = 180.0;
    private final double MIN_FIL_BRANCH_LEN = 50.0;
    private final double BRANCH_ZONE_WIDTH = RES;
    private final double T = 1.0;
    Random rand = new Random();
    private final char DELIM = '/';

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
        ip.setValue(0);
        ArrayList<Filament> filaments = new ArrayList();
        Virus virus = new Virus(RES * ip.getWidth() / 2.0, RES * ip.getHeight() / 2.0);
        for (int i = 0; i < N_FILS; i++) {
            filaments.add(createInitialFilament(virus, ip));
        }
        int filCount = 1;
        imageFolder = new File(GenUtils.openResultsDirectory((Utilities.getFolder(imageFolder,
                "Choose_Location_for_Output", true)).getAbsolutePath()
                + DELIM + "Sim_Actin_Tails" + DELIM + FIL_NOISE + "_" + +RES + "_"
                + CAP_FAC + "_" + N_FILS + "_" + P_ZONE_DEG + "_"
                + MIN_FIL_BRANCH_LEN + "_" + BRANCH_ZONE_WIDTH + "_" + virus.getBrownian()
                + "_" + virus.getMass() + "_" + virus.getcD(), "\\"));
        ProgressDialog dialog = new ProgressDialog(null, "Growing...", false, TITLE, false);
        dialog.setVisible(true);
        int virRadius = (int) Math.round(virus.getRadius() / RES);
        for (int i = 0; i < maxSteps; i++) {
            virus.brownian();
            IJ.freeMemory();
            dialog.updateProgress(i, maxSteps);
            int f1 = filCount;
            Energy netEnergy = new Energy(0.0, 0.0);
            for (int j = 0; j < f1; j++) {
                Filament current = filaments.get(j);
                if (current.grow(FIL_NOISE, virus)) {
                    int x = (int) Math.round(current.getX() / RES);
                    int y = (int) Math.round(current.getY() / RES);
                    ip.drawPixel(x, y);
                }
                double dist = Utils.calcDistance(current.getX(), current.getY(),
                        virus.getX(), virus.getY()) - virus.getRadius();
                Energy currentPE = current.calcPE(virus.getX(), virus.getY(), virus.getRadius());
                if (dist > ((rand.nextDouble() * BRANCH_ZONE_WIDTH * CAP_FAC / RES)
                        * getGrowP(virus, current.getX(), current.getY()))
                        || currentPE.getMag() > Filament.MAX_FIL_PE) {
                    filaments.remove(j);
                    j--;
                    filCount--;
                    f1--;
                    if (filCount < N_FILS) {
                        filaments.add(createInitialFilament(virus, ip));
                        filCount++;
                    }
                } else {
                    netEnergy.addEnergy(current.calcPE(virus.getX(), virus.getY(), virus.getRadius()));
                    double l = current.getLength() * current.getThickness();
                    if (l > MIN_FIL_BRANCH_LEN && rand.nextDouble() > getBranchP(dist)) {
                        double x = current.getX();
                        double y = current.getY();
                        double angle = current.getBranchAngle();
                        filaments.add(new Filament(x, y, angle));
                        current.resetLength();
                        filCount++;
                    }
                }
            }
            virus.updateVelocity(netEnergy, T);
            ImageProcessor filOut = ip.duplicate();
            ImageProcessor virOut = ip.duplicate();
            virOut.setValue(255);
            virOut.fill();
            virOut.setValue(0);
            virOut.fillOval((int) Math.round(virus.getX() / RES) - virRadius,
                    (int) Math.round(virus.getY() / RES) - virRadius,
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
        IJ.saveAs(new ImagePlus("", ip), "PNG", imageFolder + "\\Result");
        dialog.dispose();
    }

    Filament createInitialFilament(Virus virus, ImageProcessor ip) {
        double angle = Math.toRadians(getInitAngle(virus));
        double r = rand.nextDouble() * BRANCH_ZONE_WIDTH / RES + virus.getRadius();
        double x0 = virus.getX() + r * Math.cos(angle);
        double y0 = virus.getY() - r * Math.sin(angle);
        return new Filament(x0, y0, 360.0 * rand.nextDouble());
    }

    private double getInitAngle(Virus virus) {
        double avVel[] = virus.getVelAv(100);
        double angle = Utils.arcTan(avVel[0], avVel[1]) - 180.0 + P_ZONE_DEG * rand.nextGaussian();
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
        Gaussian gaussian = new Gaussian(0.0, P_ZONE_DEG);
        return gaussian.value(pos - peakPos) / gaussian.value(0.0);
    }

    private double getBranchP(double d) {
        Gaussian gaussian = new Gaussian(0.0, BRANCH_ZONE_WIDTH);
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
