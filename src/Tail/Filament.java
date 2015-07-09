package Tail;

import IAClasses.Utils;
import ij.process.FloatProcessor;
import java.util.ArrayList;
import java.util.Random;

public class Filament {

    private double x, y, branchRate;
    int length;
    ArrayList<Double> xPix = new ArrayList();
    ArrayList<Double> yPix = new ArrayList();
    Random R = new Random();
    boolean branchX = true;
    private double angle, branch = -70, branchOffset, thickness;
    Random rand = new Random();
    private double hookeK = 1.0;

    public Filament(double xc, double yc, double a0, double hgu, double thickness) {
        this.x = xc;
        this.y = yc;
        this.angle = a0;
        this.length = 0;
        this.branchRate = hgu;
        this.thickness = thickness;
        branchOffset = R.nextGaussian() * hgu * 0.2;
        if (R.nextBoolean()) {
            branchOffset *= -1;
        }
    }

    public boolean grow(double noise, Virus virus) {
        double clearance = Utils.calcDistance(x, y, virus.getX(), virus.getY());
        if (rand.nextDouble() * thickness < clearance) {
            return false;
        }
        xPix.add(x);
        yPix.add(y);
        angle += noise * R.nextGaussian();

        double xVec = thickness * Math.cos(Math.toRadians(angle));
        double yVec = thickness * Math.sin(Math.toRadians(angle));

        x += xVec;
        y += yVec;
        length++;
        return true;
    }

    Force calcForce(double xV, double yV, double r) {
        double d = (Utils.calcDistance(x, y, xV, yV) - r);
        if (Math.abs(d) > thickness) {
            return new Force(0.0, 0.0);
        } else {
            double mag = -hookeK * d;
            double xD = x - xPix.get(length - 1);
            double yD = y - yPix.get(length - 1);
            double theta = Math.PI * 2.0 * Utils.arcTan(xD, yD) / 360.0;
            return new Force(mag * Math.cos(theta), mag * Math.sin(theta));
        }
    }

    public int getLength() {
        return length;
    }

    public void resetLength() {
        Double tempX = (xPix.get(length - 1));
        Double tempY = (yPix.get(length - 1));
        xPix.clear();
        yPix.clear();
        xPix.add(tempX);
        yPix.add(tempY);
        length = 1;
    }

    public int getBranchx(double d) {
        return (int) Math.round(xPix.get(length - (int) Math.round(d) - 1));
    }

    public int getBranchy(double d) {
        return (int) Math.round(yPix.get(length - (int) Math.round(d) - 1));
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
        branchOffset = R.nextGaussian() * branchRate * 0.2;
        if (R.nextBoolean()) {
            branchOffset *= -1;
        }
    }
}
