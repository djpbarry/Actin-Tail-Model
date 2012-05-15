package Mycelium;

import ij.IJ;
import ij.ImagePlus;
import ij.ImageStack;
import ij.gui.GenericDialog;
import ij.process.ImageProcessor;
import ij.process.ImageStatistics;
import java.util.Random;

public class Mycelium {

    ImagePlus imp;
    private static double hgu = 95;
    private static double B = 0.00025;

    public void run(ImageProcessor ip) {
        int w = ip.getWidth(), h = ip.getHeight(), i, j, h0 = 0, h1, l, x0 = w / 2, y0 = h / 2, x = x0, y = y0, sw = 250, sh = 250, xlow = w / 2 - 1, xhigh = w / 2 + 1, ylow = h / 2 - 1, yhigh = h / 2 + 1;
        double angle, X = 0, apex;
        ImageStatistics is;
        Random R = new Random();

        GenericDialog gd = new GenericDialog(this.getClass().getName());
        gd.addNumericField("HGU:", hgu, 0, 6, "pixels");
        gd.addNumericField("B:", B, 5, 6, " ");
        gd.showDialog();
        if (gd.wasCanceled()) {
            return;
        }
        hgu = gd.getNextNumber();
        B = gd.getNextNumber();

        ip.setValue(0);
        ip.setLineWidth(2);
        Hypha hyphae[] = new Hypha[1000000];
        angle = R.nextDouble() * 360;

        hyphae[h0] = new Hypha(x0, y0, x0, y0, angle, ip);
        h0++;
        hyphae[h0] = new Hypha(x0, y0, x0, y0, angle + 180, ip);
        h0++;

        ImageStack stack = new ImageStack(w, h);

        while (X < B) {
            h1 = h0;
            is = imp.getStatistics();
            X = (double) is.histogram[0] / Math.pow(is.histogram[255], 1.1);
            for (j = 0; j < h1; j++) {
                hyphae[j].grow(X);
                x = hyphae[j].getX();
                y = hyphae[j].getY();
                if (x < w / 2 && x < xlow) {
                    xlow = x - 1;
                }
                if (x > w / 2 && x > xhigh) {
                    xhigh = x + 1;
                }
                if (y < h / 2 && y < ylow) {
                    ylow = y - 1;
                }
                if (y > h / 2 && y > yhigh) {
                    yhigh = y + 1;
                }
                l = hyphae[j].getLength();
                if (l > hgu + R.nextGaussian() * hgu * 0.1) {
                    apex = hgu * 0.2 + R.nextGaussian() * hgu * 0.02;
                    x = hyphae[j].getBranchx(apex);
                    y = hyphae[j].getBranchy(apex);
                    angle = hyphae[j].getBranchAngle();
                    hyphae[h0] = new Hypha(x0, y0, x, y, angle, ip);
                    hyphae[j].resetLength();
                    h0++;
                }
            }
            stack.addSlice("" + X, ip.duplicate());
            IJ.showProgress(X / 0.05);
        }
        ImagePlus anim = new ImagePlus("Plot", stack);
        anim.show();
    }

    public class Hypha {

        private double x0, y0, x, y, angle, branch = -90;
        int dim, length;
        double xpix[] = new double[(int) (hgu * 1.2)];
        double ypix[] = new double[(int) (hgu * 1.2)];
        ImageProcessor plot;
        Random R = new Random();

        public Hypha(double x1, double y1, double xc, double yc, double a0, ImageProcessor ip) {
            x0 = x1;
            y0 = y1;
            x = xc;
            y = yc;
            angle = a0;
            plot = ip;
            length = 0;
            dim = Math.min(ip.getWidth(), ip.getHeight());
        }

        public void grow(double biomass) {
            xpix[length] = x;
            ypix[length] = y;

            double m = (y - y0) / (x0 - x);
            double t = Math.toDegrees(Math.atan(m));
            if (x < x0) {
                t += 180;
            }

            if (x == x0 && y <= y0) {
                t = 90;
            }
            if (x == x0 && y >= y0) {
                t = 270;
            }

            if (t < 0) {
                t += 360;
            }

            angle += biomass * (t - angle) + 10 * R.nextGaussian();

            x += Math.cos(Math.toRadians(angle));
            y -= Math.sin(Math.toRadians(angle));
            plot.drawLine((int) xpix[length], (int) ypix[length], (int) x, (int) y);
            length++;
        }

        public int getLength() {
            return length;
        }

        public void resetLength() {
            length = 0;
        }

        public double getBranchAngle() {
            branch += -2 * branch;
            return (branch + angle);
        }

        public int getBranchx(double d) {
            return (int) xpix[length - (int) d];
        }

        public int getBranchy(double d) {
            return (int) ypix[length - (int) d];
        }

        public int getX() {
            return (int) x;
        }

        public int getY() {
            return (int) y;
        }
    }
}
