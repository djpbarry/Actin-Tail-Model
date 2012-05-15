package Mycelium;

import ij.IJ;
import ij.process.ImageProcessor;

public class Hypha {

    private double x0, y0, x, y, angle, branch = -90;
    private int length;
    ImageProcessor plot;

    public Hypha(double x1, double y1, double xc, double yc, double a0, ImageProcessor ip) {
        x0 = x1;
        y0 = y1;
        x = xc;
        y = yc;
        angle = a0;
        plot = ip;
        length = 0;
    }

    public void grow(double biomass) {
        double oldx = x;
        double oldy = y;

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

        IJ.write("" + m + " " + t);

        x += Math.cos(Math.toRadians(angle));
        y -= Math.sin(Math.toRadians(angle));
        plot.drawLine((int) oldx, (int) oldy, (int) x, (int) y);
        length++;
    }

    public int getLength() {
        return length;
    }

    public double getBranchAngle() {
        branch += -2 * branch;
        return (branch + angle);
    }

    public int getBranchx(int d) {
        return (int) Math.round(x -= d * Math.cos(Math.toRadians(angle)));
    }

    public int getBranchy(int d) {
        return (int) Math.round(y += d * Math.sin(Math.toRadians(angle)));
    }
}
