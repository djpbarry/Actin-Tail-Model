package Tail;

import IAClasses.ProgressDialog;
import IAClasses.Utils;
import UtilClasses.GenUtils;
import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.process.ByteProcessor;
import ij.process.ColorProcessor;
import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Date;
import java.util.LinkedList;
import java.util.Random;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.stat.descriptive.summary.Sum;

public class TailGrower {

    private final double FIL_NOISE = 5.0;
    public File imageFolder, parentFolder = new File("c:/users/barry05/desktop/Sim_Actin_Tails");
    private static int maxSteps = 4000, width = 1000, height = 1000;
    private final String TITLE = "";
    private static boolean showAllImages = true;
    private final double RES = 7; // One pixel equals 7 nm, approximate width of actin filament - http://www.ncbi.nlm.nih.gov/books/NBK9908/
    private final double CAP_FAC = 150.0;
    private final int N_FILS = 20;
    private final int NPFS = 20;
    private final double SIGMA = RES / 2.0;
//    private final double P_ZONE_DEG;
    private final int P_ZONE_STEP = 5, P_ZONE_MAX = 130;
    private final double MIN_FIL_BRANCH_LEN = 50.0;
    private final double BRANCH_ZONE_WIDTH = RES;
    private final double T = 1.0 / 250.0;
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
        parentFolder = new File(GenUtils.openResultsDirectory(parentFolder.getAbsolutePath()
                + DELIM + ((new Date()).toString()).replace(':', '-'), "/"));
        try {
            velFile = new File(parentFolder + "/velocities.csv");
            velstream = new PrintWriter(new FileOutputStream(velFile));

        } catch (FileNotFoundException e) {
            System.out.println(e.toString());
            return;
        }
        velstream.println("FIL_NOISE," + FIL_NOISE + ",RES," + RES + ",CAP_FAC,"
                + CAP_FAC + ",N_FILS," + N_FILS + ",MIN_FIL_BRANCH_LEN,"
                + MIN_FIL_BRANCH_LEN + ",BRANCH_ZONE_WIDTH," + BRANCH_ZONE_WIDTH);
        for (int pzone = P_ZONE_STEP; pzone <= P_ZONE_MAX; pzone += P_ZONE_STEP) {
            ColorProcessor bp = new ColorProcessor(width, height);
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

    public void grow(ColorProcessor ip, int maxSteps, double[] vels, double pZone) {
        ByteProcessor result = new ByteProcessor(ip.getWidth(), ip.getHeight());
        result.setColor(Color.white);
        result.fill();
        result.setColor(Color.black);
        DecimalFormat numFormat = new DecimalFormat("000");
        ArrayList<Filament> filaments = new ArrayList();
        double xv0 = RES * ip.getWidth() / 2.0, yv0 = RES * ip.getHeight() / 2.0;
        double dists[] = new double[maxSteps];
//        ArrayList<Color> colours = new ArrayList();
//        double vels[] = new double[maxSteps];
        Virus virus = new Virus(xv0, yv0);
        for (int i = 0; i < N_FILS; i++) {
            filaments.add(createInitialFilament(virus, pZone, 0));
//            colours.add(new Color(rand.nextInt(256), rand.nextInt(256), rand.nextInt(256)));
        }
//        imageFolder = new File(GenUtils.openResultsDirectory((Utilities.getFolder(imageFolder,
//                "Choose_Location_for_Output", true)).getAbsolutePath()
//                + DELIM + "Sim_Actin_Tails" + DELIM + "Output", "\\"));
        imageFolder = new File(GenUtils.openResultsDirectory(parentFolder.getAbsolutePath()
                + DELIM + pZone, "\\"));
        File traj, animData;
        PrintWriter trajStream, animStream;
        try {
            traj = new File(imageFolder + "/trajectory.csv");
            trajStream = new PrintWriter(new FileOutputStream(traj));
            animData = new File(imageFolder + "/anim_data.txt");
            animStream = new PrintWriter(new FileOutputStream(animData));
        } catch (FileNotFoundException e) {
            System.out.println(e.toString());
            return;
        }
        trajStream.println("FIL_NOISE," + FIL_NOISE + ",RES," + RES + ",CAP_FAC,"
                + CAP_FAC + ",N_FILS," + N_FILS + ",P_ZONE_DEG," + pZone + ",MIN_FIL_BRANCH_LEN,"
                + MIN_FIL_BRANCH_LEN + ",BRANCH_ZONE_WIDTH," + BRANCH_ZONE_WIDTH + ",VIRUS_BROWNIAN," + virus.getBrownian()
                + ",VIRUS_MASS," + virus.getMass() + ",VIRUS_CD," + virus.getcD());
        trajStream.println("t,x,y,xVel,yVel");
        animStream.println("FRAMES " + maxSteps);
        animStream.println("WIDTH " + ip.getWidth());
        animStream.println("HEIGHT " + ip.getHeight());
        ProgressDialog dialog = new ProgressDialog(null, "Growing...", false, TITLE, false);
        dialog.setVisible(true);
        int virRadius = (int) Math.round(virus.getRadius() / RES);
        for (int i = 0; i < maxSteps; i++) {
            int uncapped = 0;
            ip = new ColorProcessor(ip.getWidth(), ip.getHeight());
            ip.setColor(Color.white);
            ip.fill();
            trajStream.println(i + "," + virus.getX() + "," + virus.getY() + "," + virus.getxVel() + "," + virus.getyVel());
            virus.brownian();
            animStream.print(virus.getX() + " " + virus.getY() + " ");
            IJ.freeMemory();
            dialog.updateProgress(i, maxSteps);
            int f1 = filaments.size();
            Energy netEnergy = new Energy(0.0, 0.0, 0.0);
            for (int j = 0; j < f1; j++) {
                Filament current = filaments.get(j);
                if (!current.isCapped()) {
                    uncapped++;
                    Monomer end = current.getMons().getLast();
                    if (current.grow(FIL_NOISE, virus, filaments, j, i)) {
                        end = current.getMons().getLast();
                        animStream.print(end.getX() + " " + end.getY() + " ");
                    }
                    double dist = Utils.calcDistance(end.getX(), end.getY(),
                            virus.getX(), virus.getY()) - virus.getRadius();
                    Energy currentPE = current.calcRPE(virus.getX(), virus.getY(), virus.getRadius());
                    if (dist > ((rand.nextDouble() * BRANCH_ZONE_WIDTH * CAP_FAC / RES)
                            * getGrowP(virus, end.getX(), end.getY(), pZone, SIGMA, NPFS))
                            || currentPE.getMag() > Filament.MAX_FIL_PE) {
                        current.setCapped(true);
                    } else {
                        netEnergy.addEnergy(current.calcRPE(virus.getX(), virus.getY(), virus.getRadius()));
                        double l = current.getLength() * current.getThickness();
                        if (l > MIN_FIL_BRANCH_LEN && rand.nextDouble() > getBranchP(dist)) {
                            double x = end.getX();
                            double y = end.getY();
                            double angle = current.getBranchAngle();
                            filaments.add(new Filament(x, y, angle, i, T));
                            uncapped++;
//                        colours.add(new Color(rand.nextInt(256), rand.nextInt(256), rand.nextInt(256)));
                            current.resetLength();
                        }
                    }
                }
                current.updateMonomerStates(i);
                LinkedList<Monomer> mons = current.getMons();
                if (mons.size() < 1) {
                    filaments.remove(j);
//                    colours.remove(j);
                    j--;
                    f1--;
                } else {
//                    ip.setColor(labels.get(j));
                    for (Monomer m : mons) {
                        int x = (int) Math.round(m.getX() / RES);
                        int y = (int) Math.round(m.getY() / RES);
                        if (m.getState() == Monomer.ADP) {
                            ip.setColor(Color.blue);
                        } else if (m.getState() == Monomer.ATPPI) {
                            ip.setColor(Color.green);
                        } else {
                            ip.setColor(Color.red);
                        }
                        ip.drawPixel(x, y);
                        result.drawPixel(x, y);
                    }
                }
            }
            virus.updateVelocity(netEnergy, T);
            if (showAllImages) {
//                ImageProcessor filOut = ip.duplicate();
                ByteProcessor virOut = new ByteProcessor(ip.getWidth(), ip.getHeight());
                virOut.setValue(255);
                virOut.fill();
                virOut.setValue(0);
                double theta = Math.toRadians(virus.getTheta());
                virOut.drawOval((int) Math.round(virus.getX() / RES - virRadius),
                        (int) Math.round(virus.getY() / RES - virRadius), 2 * virRadius, 2 * virRadius);
                virOut.drawLine((int) Math.round(virus.getX() / RES + virRadius * Math.cos(theta)),
                        (int) Math.round(virus.getY() / RES - virRadius * Math.sin(theta)),
                        (int) Math.round(virus.getX() / RES - virRadius * Math.cos(theta)),
                        (int) Math.round(virus.getY() / RES + virRadius * Math.sin(theta)));
//                FloatProcessor growthZone = new FloatProcessor(ip.getWidth(), ip.getHeight());
//                growthZone.setValue(0);
//                growthZone.fill();
//                double x0 = virus.getX();
//                double y0 = virus.getY();
//                double r = 2.0 * virus.getRadius();
//                for (double j = y0 - r; j <= y0 + r; j += RES) {
//                    for (double k = x0 - r; k <= x0 + r; k += RES) {
//                        growthZone.putPixelValue((int) Math.round(k / RES),
//                                (int) Math.round(j / RES), getGrowP(virus, k, j, pZone, SIGMA, NPFS));
//                    }
//                }
                IJ.saveAs(new ImagePlus("", ip.duplicate()), "PNG", imageFolder + "\\Filaments_Step" + numFormat.format(i));
                IJ.saveAs(new ImagePlus("", virOut), "PNG", imageFolder + "\\Virus_Step" + numFormat.format(i));
//                IJ.saveAs(new ImagePlus("", growthZone), "tif", imageFolder + "\\GrowthZone_Step" + numFormat.format(i));
            }
            vels[i] = Math.sqrt(Math.pow(virus.getxVel(), 2.0) + Math.pow(virus.getyVel(), 2.0));
            dists[i] = vels[i] * T;
            animStream.println();
            if (uncapped < N_FILS) {
                //TODO
                filaments.add(createInitialFilament(virus, pZone, i));
//                    colours.add(new Color(rand.nextInt(256), rand.nextInt(256), rand.nextInt(256)));
            }
        }
        trajStream.close();
        animStream.close();
        IJ.saveAs(new ImagePlus("", result), "PNG", imageFolder + "\\Result");
        dialog.dispose();
        Sum sum = new Sum();
        System.out.println(pZone + ","
                + (Utils.calcDistance(xv0, yv0, virus.getX(), virus.getY()) / sum.evaluate(dists, 0, maxSteps)));
//        printStats(vels);
    }

    Filament createInitialFilament(Virus virus, double pZone, int time) {
        double angle = Math.toRadians(getInitAngle(virus, pZone));
        double r = rand.nextDouble() * BRANCH_ZONE_WIDTH / RES + virus.getRadius();
        double x0 = virus.getX() + r * Math.cos(angle);
        double y0 = virus.getY() - r * Math.sin(angle);
        return new Filament(x0, y0, 360.0 * rand.nextDouble(), time, T);
    }

    private double getInitAngle(Virus virus, double P_ZONE_DEG) {
        double avVel[] = virus.getVelAv();
        double angle = Utils.arcTan(avVel[0], avVel[1]) - 180.0 + P_ZONE_DEG * rand.nextGaussian();
        if (angle < 0.0) {
            angle += 360.0;
        }
        return angle;
    }

    private double getGrowP(Virus virus, double x, double y, double pZone, double sigma, int sites) {
        double peakPos = virus.getTheta();
        double pos = Utils.arcTan(x - virus.getX(), y - virus.getY());
        double diff = pos - peakPos;
        if (diff < 0.0) {
            diff += 360.0;
        }
        if (diff > 180.0) {
            diff = 360.0 - diff;
        }
//        Gaussian gaussian = new Gaussian(0.0, pZone);
        return getNPFP(diff, sigma, pZone, sites) / getNPFP(0.0, sigma, pZone, sites);
    }

    double getNPFP(double x, double sigma, double pZone, int sites) {
        int S = (sites - 1) / 2;
        double thetaStep = pZone / S;
        double sum = 0.0;
        for (int i = 0; i <= S; i++) {
            double thisX = 0.0 + i * thetaStep;
            sum += (new Gaussian(thisX, sigma)).value(x);
        }
        return sum;
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
