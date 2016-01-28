package Tail;

import IAClasses.Utils;
import java.util.ArrayList;
import java.util.LinkedList;
import java.util.Random;

public class Filament {

    int length;
    LinkedList<Monomer> mons = new LinkedList();
    Random rand = new Random();
    boolean branchX = true;
    private double angle, branch = -70;
    private double hookeFil = 10.0, hookeBond = 1.0;
    public final static double MAX_FIL_PE = 1000.0, MAX_BOND_PE = 5.0;

    public Filament(double xc, double yc, double a0) {
        mons.add(new Monomer(xc, yc));
        this.angle = a0;
        this.length = 1;
        if (rand.nextBoolean()) {
            branch *= -1;
        }
    }

    public boolean grow(double noise, Virus virus, ArrayList<Filament> filaments, int index) {
        Monomer m = mons.getLast();
        double x = m.getX();
        double y = m.getY();
        double minClearance = rand.nextDouble() * Monomer.DIAMETER;
        double virClearance = Utils.calcDistance(x, y, virus.getX(), virus.getY()) - virus.getRadius();
        double filClearance = getMinFilClearance(filaments, index);
        if (virClearance < minClearance || filClearance < minClearance) {
            return false;
        }
        angle += noise * rand.nextGaussian();

        double xVec = Monomer.DIAMETER * Math.cos(Math.toRadians(angle));
        double yVec = -Monomer.DIAMETER * Math.sin(Math.toRadians(angle));
        mons.add(new Monomer(x + xVec, y + yVec));
        length++;
        return true;
    }

    double getMinFilClearance(ArrayList<Filament> filaments, int index) {
        double minDist = Double.MAX_VALUE;
        for (int i = 0; i < index; i++) {
            Filament current = filaments.get(i);
            Monomer m1 = mons.getLast();
            Monomer m2 = current.getMons().getLast();
            double dist = Utils.calcDistance(m2.getX(), m2.getY(), m1.getX(), m2.getY());
            if (dist < minDist) {
                minDist = dist;
            }
        }
        return minDist;
    }

    Energy calcRPE(double xV, double yV, double r) {
        Monomer m = mons.getLast();
        double x = m.getX();
        double y = m.getY();
        double d = Utils.calcDistance(x, y, xV, yV) - r;
        double E;
        if (d < 0.0) {
            E = 0.5 * hookeFil * Math.pow(d, 2.0);
        } else if (d > 0.0 && d < Monomer.DIAMETER) {
            E = -0.5 * hookeBond * Math.pow(d, 2.0);
        } else {
            return new Energy(0.0, 0.0, 0.0);
        }
        double xD = xV - x;
        double yD = yV - y;
        double gamma = Math.toRadians(Utils.arcTan(-xD, -yD));
        double dir[] = {Monomer.DIAMETER * Math.cos(Math.toRadians(angle)), -Monomer.DIAMETER * Math.sin(Math.toRadians(angle))};
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
        Monomer temp = (mons.get(length - 1));
        mons.clear();
        mons.add(temp);
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

    public LinkedList<Monomer> getMons() {
        return mons;
    }

    public double getThickness() {
        return Monomer.DIAMETER;
    }
}
