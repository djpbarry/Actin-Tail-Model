package Tail;

import IAClasses.Utils;
import java.util.ArrayList;
import java.util.Random;

public class Filament {

    private double x, y;
    int length;
    ArrayList<Double> xPix = new ArrayList();
    ArrayList<Double> yPix = new ArrayList();
    Random rand = new Random();
    boolean branchX = true;
    private double angle, branch = -70, thickness = 7;
    private double hookeFil = 10.0, hookeBond = 1.0;

    public Filament(double xc, double yc, double a0) {
        this.x = xc;
        this.y = yc;
        this.angle = a0;
        this.length = 0;
        if (rand.nextBoolean()) {
            branch *= -1;
        }
    }

    public boolean grow(double noise, Virus virus) {
        double clearance = Utils.calcDistance(x, y, virus.getX(), virus.getY()) - virus.getRadius();
        if (clearance < rand.nextDouble() * thickness) {
            return false;
        }
        xPix.add(x);
        yPix.add(y);
        angle += noise * rand.nextGaussian();

        double xVec = thickness * Math.cos(Math.toRadians(angle));
        double yVec = -thickness * Math.sin(Math.toRadians(angle));

        x += xVec;
        y += yVec;
        length++;
        return true;
    }

    Energy calcPE(double xV, double yV, double r) {
        double d = Utils.calcDistance(x, y, xV, yV) - r;
        double mag = 0.0;
        if (d < 0.0) {
            mag = 0.5 * hookeFil * Math.pow(d, 2.0);
        } else if (d > 0.0 && d < thickness) {
            mag = -0.5 * hookeBond * Math.pow(d, 2.0);
        }
        double xD = xV - x;
        double yD = yV - y;
        double theta = Math.toRadians(Utils.arcTan(xD, yD));
        return new Energy(mag * Math.cos(theta), -mag * Math.sin(theta));
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

    public double getBranchAngle() {
        branch += -2 * branch;
        return (branch + angle);
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getThickness() {
        return thickness;
    }
}
