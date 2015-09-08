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
    public final static double MAX_FIL_PE = 100.0, MAX_BOND_PE = 5.0;

    public Filament(double xc, double yc, double a0) {
        this.x = xc;
        this.y = yc;
        xPix.add(x);
        yPix.add(y);
        this.angle = a0;
        this.length = 1;
        if (rand.nextBoolean()) {
            branch *= -1;
        }
    }

    public boolean grow(double noise, Virus virus) {
        double clearance = Utils.calcDistance(x, y, virus.getX(), virus.getY()) - virus.getRadius();
        if (clearance < rand.nextDouble() * thickness) {
            return false;
        }
        angle += noise * rand.nextGaussian();

        double xVec = thickness * Math.cos(Math.toRadians(angle));
        double yVec = -thickness * Math.sin(Math.toRadians(angle));

        x += xVec;
        y += yVec;
        xPix.add(x);
        yPix.add(y);
        length++;
        return true;
    }

    Energy calcRPE(double xV, double yV, double r) {
        double d = Utils.calcDistance(x, y, xV, yV) - r;
        double E;
        if (d < 0.0) {
            E = 0.5 * hookeFil * Math.pow(d, 2.0);
        } else if (d > 0.0 && d < thickness) {
            E = -0.5 * hookeBond * Math.pow(d, 2.0);
        } else {
            return new Energy(0.0, 0.0, 0.0);
        }
        double xD = xV - x;
        double yD = yV - y;
        double gamma = Math.toRadians(Utils.arcTan(-xD, -yD));
        double dir[] = {thickness * Math.cos(Math.toRadians(angle)), -thickness * Math.sin(Math.toRadians(angle))};
        double yFt2 = (x - xV) * Math.sin(gamma) + (y - yV) * Math.cos(gamma) + yV;
        double yFt1 = (x - dir[0] - xV) * Math.sin(gamma) + (y - dir[0] - yV) * Math.cos(gamma) + yV;
        double theta = Math.PI / 2.0 - Utils.angleBetweenTwoLines(new double[]{Math.cos(angle), -Math.sin(angle)}, new double[]{xD, yD});
        double Er = E * Math.sin(theta);
        double Et = E * Math.cos(theta);
        if (yFt2 > yFt1) {
            Et *= -1.0;
        }
        double phi = Math.toRadians(Utils.arcTan(xD, yD));
        return new Energy(Er * Math.cos(phi), -Er * Math.sin(phi), Et);
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
        double a = branch + angle;
        if (a < 0.0) {
            a += 360.0;
        } else if (a >= 360.0) {
            a -= 360.0;
        }
        return a;
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
