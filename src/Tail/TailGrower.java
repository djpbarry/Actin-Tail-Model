package Tail;

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
import java.util.Date;
import java.util.Random;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.stat.descriptive.summary.Sum;

public class TailGrower {

    private final double FIL_NOISE = 5.0;
    public File imageFolder, parentFolder = new File("c:/users/barry05/desktop/Sim_Actin_Tails");
    private static int maxSteps = 1000, width = 1000, height = 1000;
    private final String TITLE = "";
    private static boolean showAllImages = true;
    private final double RES = 7; // One pixel equals 7 nm, approximate width of actin filament - http://www.ncbi.nlm.nih.gov/books/NBK9908/
    private final double CAP_FAC = 150.0;
    private final int N_FILS = 20;
//    private final double P_ZONE_DEG;
    private final int P_ZONE_STEP = 5, P_ZONE_MAX = 180;
    private final double MIN_FIL_BRANCH_LEN = 50.0;
    private final double BRANCH_ZONE_WIDTH = RES;
    private final double T = 1.0;
    Random rand = new Random();
    private final char DELIM = '/';

    public static void main(String args[]) {
//        for (int pzone = 5; pzone <= 180; pzone += 5) {
        TailGrower grower = new TailGrower();
//            if (!grower.initialise()) {
//                System.exit(-1);
//            }
        grower.run();
//        }
        System.exit(0);
    }

    boolean initialise() {
        if (!showDialog(TITLE)) {
            return false;
        }
        return true;
    }

    void run() {
        double vels[][] = new double[P_ZONE_MAX / P_ZONE_STEP + 1][maxSteps];
        File velFile;
        PrintWriter velstream;
        parentFolder = new File(GenUtils.openResultsDirectory(parentFolder.getAbsolutePath(), "/"));
        try {
            velFile = new File(parentFolder + "/velocities.csv");
            velstream = new PrintWriter(new FileOutputStream(velFile));
            velstream.println("FIL_NOISE," + FIL_NOISE + ",RES," + RES + ",CAP_FAC,"
                    + CAP_FAC + ",N_FILS," + N_FILS + ",MIN_FIL_BRANCH_LEN,"
                    + MIN_FIL_BRANCH_LEN + ",BRANCH_ZONE_WIDTH," + BRANCH_ZONE_WIDTH);
        } catch (FileNotFoundException e) {
            System.out.println(e.toString());
            return;
        }
        for (int pzone = P_ZONE_STEP; pzone <= P_ZONE_MAX; pzone += P_ZONE_STEP) {
            ByteProcessor bp = new ByteProcessor(width, height);
            bp.setColor(Color.white);
            bp.fill();
            grow(bp, maxSteps, vels[pzone / P_ZONE_STEP - 1], pzone);
             velstream.print(pzone + ",");
        }
        velstream.println();
        for (int i = 0; i < maxSteps; i++) {
            for (int pzone = P_ZONE_STEP; pzone <= P_ZONE_MAX; pzone += P_ZONE_STEP) {
                velstream.print(vels[pzone / P_ZONE_STEP - 1][i] + ",");
            }
            velstream.println();
        }
        velstream.close();
    }

    public TailGrower() {
    }

    public void grow(ByteProcessor ip, int maxSteps, double[] vels, double pZone) {
        Date d = new Date();
        String date = d.toString();
        date = date.replace(':', '-');
        DecimalFormat numFormat = new DecimalFormat("000");
        ip.setValue(0);
        ArrayList<Filament> filaments = new ArrayList();
        double x0 = RES * ip.getWidth() / 2.0, y0 = RES * ip.getHeight() / 2.0;
        double dists[] = new double[maxSteps];
//        ArrayList<Color> colours = new ArrayList();
//        double vels[] = new double[maxSteps];
        Virus virus = new Virus(x0, y0);
        for (int i = 0; i < N_FILS; i++) {
            filaments.add(createInitialFilament(virus, ip, pZone));
//            colours.add(new Color(rand.nextInt(256), rand.nextInt(256), rand.nextInt(256)));
        }
        int filCount = 1;
//        imageFolder = new File(GenUtils.openResultsDirectory((Utilities.getFolder(imageFolder,
//                "Choose_Location_for_Output", true)).getAbsolutePath()
//                + DELIM + "Sim_Actin_Tails" + DELIM + "Output", "\\"));
        imageFolder = new File(GenUtils.openResultsDirectory(parentFolder.getAbsolutePath()
                + DELIM + date + DELIM + pZone, "\\"));
        File traj;
        PrintWriter trajStream;
        try {
            traj = new File(imageFolder + "/trajectory.csv");
            trajStream = new PrintWriter(new FileOutputStream(traj));
            trajStream.println("FIL_NOISE," + FIL_NOISE + ",RES," + RES + ",CAP_FAC,"
                    + CAP_FAC + ",N_FILS," + N_FILS + ",P_ZONE_DEG," + pZone + ",MIN_FIL_BRANCH_LEN,"
                    + MIN_FIL_BRANCH_LEN + ",BRANCH_ZONE_WIDTH," + BRANCH_ZONE_WIDTH + ",VIRUS_BROWNIAN," + virus.getBrownian()
                    + ",VIRUS_MASS," + virus.getMass() + ",VIRUS_CD," + virus.getcD());
            trajStream.println("t,x,y,xVel,yVel");
        } catch (FileNotFoundException e) {
            System.out.println(e.toString());
            return;
        }
        ProgressDialog dialog = new ProgressDialog(null, "Growing...", false, TITLE, false);
        dialog.setVisible(true);
        int virRadius = (int) Math.round(virus.getRadius() / RES);
        for (int i = 0; i < maxSteps; i++) {
            trajStream.println(i + "," + virus.getX() + "," + virus.getY() + "," + virus.getxVel() + "," + virus.getyVel());
            virus.brownian();
            IJ.freeMemory();
            dialog.updateProgress(i, maxSteps);
            int f1 = filCount;
            Energy netEnergy = new Energy(0.0, 0.0, 0.0);
            for (int j = 0; j < f1; j++) {
                Filament current = filaments.get(j);
                if (current.grow(FIL_NOISE, virus, filaments, j)) {
                    int x = (int) Math.round(current.getX() / RES);
                    int y = (int) Math.round(current.getY() / RES);
                    ip.drawPixel(x, y);
                }
                double dist = Utils.calcDistance(current.getX(), current.getY(),
                        virus.getX(), virus.getY()) - virus.getRadius();
                Energy currentPE = current.calcRPE(virus.getX(), virus.getY(), virus.getRadius());
                if (dist > ((rand.nextDouble() * BRANCH_ZONE_WIDTH * CAP_FAC / RES)
                        * getGrowP(virus, current.getX(), current.getY(), pZone))
                        || currentPE.getMag() > Filament.MAX_FIL_PE) {
                    filaments.remove(j);
//                    colours.remove(j);
                    j--;
                    filCount--;
                    f1--;
                    if (filCount < N_FILS) {
                        //TODO
                        filaments.add(createInitialFilament(virus, ip, pZone));
//                        colours.add(new Color(rand.nextInt(256), rand.nextInt(256), rand.nextInt(256)));
                        filCount++;
                    }
                } else {
                    netEnergy.addEnergy(current.calcRPE(virus.getX(), virus.getY(), virus.getRadius()));
                    double l = current.getLength() * current.getThickness();
                    if (l > MIN_FIL_BRANCH_LEN && rand.nextDouble() > getBranchP(dist)) {
                        double x = current.getX();
                        double y = current.getY();
                        double angle = current.getBranchAngle();
                        filaments.add(new Filament(x, y, angle));
//                        colours.add(new Color(rand.nextInt(256), rand.nextInt(256), rand.nextInt(256)));
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
//            double theta = Math.toRadians(virus.getTheta());
            virOut.fillOval((int) Math.round(virus.getX() / RES - virRadius),
                    (int) Math.round(virus.getY() / RES - virRadius), 2 * virRadius, 2 * virRadius);
//            virOut.drawLine((int) Math.round(virus.getX() / RES + virRadius * Math.cos(theta)),
//                    (int) Math.round(virus.getY() / RES - virRadius * Math.sin(theta)),
//                    (int) Math.round(virus.getX() / RES - virRadius * Math.cos(theta)),
//                    (int) Math.round(virus.getY() / RES + virRadius * Math.sin(theta)));
//            ByteProcessor growthZone = new ByteProcessor(ip.getWidth(), ip.getHeight());
//            growthZone.setValue(0);
//            growthZone.fill();
//            double x0 = virus.getX();
//            double y0 = virus.getY();
//            double r = 2.0 * virus.getRadius();
//            for (double j = y0 - r; j <= y0 + r; j += RES) {
//                for (double k = x0 - r; k <= x0 + r; k += RES) {
//                    growthZone.putPixelValue((int) Math.round(k / RES),
//                            (int) Math.round(j / RES), (int) Math.round(255.0 * getGrowP(virus, k, j)));
//                }
//            }
            if (showAllImages) {
                IJ.saveAs(new ImagePlus("", filOut), "PNG", imageFolder + "\\Filaments_Step" + numFormat.format(i));
                IJ.saveAs(new ImagePlus("", virOut), "PNG", imageFolder + "\\Virus_Step" + numFormat.format(i));
//                IJ.saveAs(new ImagePlus("", growthZone), "PNG", imageFolder + "\\GrowthZone_Step" + numFormat.format(i));
            }
            vels[i] = Math.sqrt(Math.pow(virus.getxVel(), 2.0) + Math.pow(virus.getyVel(), 2.0));
            dists[i] = vels[i] * T;
        }
        trajStream.close();
        IJ.saveAs(new ImagePlus("", ip), "PNG", imageFolder + "\\Result");
        dialog.dispose();
        Sum sum = new Sum();
        System.out.println(pZone + ","
                + (Utils.calcDistance(x0, y0, virus.getX(), virus.getY()) / sum.evaluate(dists, 0, maxSteps)));
//        printStats(vels);
    }

    Filament createInitialFilament(Virus virus, ImageProcessor ip, double pZone) {
        double angle = Math.toRadians(getInitAngle(virus, pZone));
        double r = rand.nextDouble() * BRANCH_ZONE_WIDTH / RES + virus.getRadius();
        double x0 = virus.getX() + r * Math.cos(angle);
        double y0 = virus.getY() - r * Math.sin(angle);
        return new Filament(x0, y0, 360.0 * rand.nextDouble());
    }

    private double getInitAngle(Virus virus, double P_ZONE_DEG) {
        double avVel[] = virus.getVelAv();
        double angle = Utils.arcTan(avVel[0], avVel[1]) - 180.0 + P_ZONE_DEG * rand.nextGaussian();
        if (angle < 0.0) {
            angle += 360.0;
        }
        return angle;
    }

    private double getGrowP(Virus virus, double x, double y, double pZone) {
        double peakPos = virus.getTheta();
        double pos = Utils.arcTan(x - virus.getX(), y - virus.getY());
        double diff = pos - peakPos;
        if (diff < 0.0) {
            diff += 360.0;
        }
        if (diff > 180.0) {
            diff = 360.0 - diff;
        }
        Gaussian gaussian = new Gaussian(0.0, pZone);
        return gaussian.value(diff) / gaussian.value(0.0);
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

//    void printStats(double[] vels) {
//        Mean mean = new Mean();
//        StandardDeviation std = new StandardDeviation();
//        System.out.println(P_ZONE_DEG + "," + mean.evaluate(vels, 0, vels.length) + "," + std.evaluate(vels, 0, vels.length));
//    }
}
