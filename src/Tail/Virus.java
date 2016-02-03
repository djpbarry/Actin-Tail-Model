/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Tail;

import java.util.ArrayList;
import java.util.Random;

/**
 *
 * @author David Barry <david.barry at crick.ac.uk>
 */
public class Virus {

    private double x, y, theta, xVel, yVel, wVel;
//    private final double mass = 10e-15; // Mass of single vaccinia virion is ~10 fg - dx.doi.org/10.1016/j.snb.2005.08.047
    private final double mass = 0.4;
    private final double rho = 1.0/150.0;
    private final double radius = 150.0; // Width of virion is 250 - 350nm
    private final double area = Math.PI * radius * radius;
    private final double cD = 0.0001;
    private final double mu = 40.0;
    private final double inertia = mass * Math.pow(radius, 3.0) / 3;
    private Random r = new Random();
    private ArrayList<double[]> vels;
    private double brownian = 10.0;
    private final int STEP_SIZE = 100;

    public Virus(double x, double y) {
        this.x = x;
        this.y = y;
        xVel = 0.0;
        yVel = 0.0;
        wVel = 0.0;
        vels = new ArrayList();
        vels.add(new double[]{xVel, yVel});
        theta = r.nextDouble() * 360.0;
    }

    public void updateVelocity(Energy Es, double t) {
        double xInc = 2.0 * Es.getxE() / mass;
        if (Es.getxE() > 0.0) {
            xVel += Math.sqrt(xInc);
        } else {
            xVel += -Math.sqrt(-xInc);
        }
        double xDrag = t * calcDragForce(xVel) / mass;
        if (Math.abs(xDrag) > Math.abs(xVel)) {
            xVel = 0.0;
        } else {
            xVel -= t * calcDragForce(xVel) / mass;
        }
        double yInc = 2.0 * Es.getyE() / mass;
        if (Es.getyE() > 0.0) {
            yVel += Math.sqrt(yInc);
        } else {
            yVel += -Math.sqrt(-yInc);
        }
        double yDrag = t * calcDragForce(yVel) / mass;
        if (Math.abs(yDrag) > Math.abs(yVel)) {
            yVel = 0.0;
        } else {
            yVel -= t * calcDragForce(yVel) / mass;
        }
        double wInc;
        if (Es.getwE() > 0.0) {
            wInc = Math.sqrt(2.0 * Es.getwE() / inertia);
        } else {
            wInc = -Math.sqrt(-2.0 * Es.getwE() / inertia);
        }
        wVel += wInc;
        double wDrag = t * calcRotFricForce(wVel) / mass;
        if (Math.abs(wDrag) > Math.abs(wVel)) {
            wVel = 0.0;
        } else {
            wVel -= wDrag;
        }
        vels.add(new double[]{xVel, yVel});
        x += xVel * t;
        y += yVel * t;
        updateTheta(Math.toDegrees(wVel * t));
    }

    public void brownian() {
        x += -brownian / 2.0 + brownian * r.nextDouble();
        y += -brownian / 2.0 + brownian * r.nextDouble();
        updateTheta(-brownian / 2.0 + brownian * r.nextDouble());
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getTheta() {
        return theta;
    }

    public double getRadius() {
        return radius;
    }

    double calcDragForce(double vel) {
        double vInc = 0.5 * rho * area * cD * vel * vel;
        if (vel < 0.0) {
            vInc *= -1.0;
        }
        return vInc;
    }

    double calcRotFricForce(double vel) {
        return mu * vel;
    }

    public double[] getVelAv() {
        double xVA = 0.0;
        double yVA = 0.0;
        int s = vels.size();
        for (int i = 0; i < STEP_SIZE && i < s; i++) {
            xVA += vels.get(s - i - 1)[0];
            yVA += vels.get(s - i - 1)[1];
        }
        return new double[]{xVA / STEP_SIZE, yVA / STEP_SIZE};
    }

    public double getMass() {
        return mass;
    }

    public double getcD() {
        return cD;
    }

    public double getBrownian() {
        return brownian;
    }

    void updateTheta(double inc) {
        theta += inc;
        if (theta < 0.0) {
            theta += 360.0;
        } else if (theta >= 360.0) {
            theta -= 360;
        }
    }

    public double getxVel() {
        return xVel;
    }

    public double getyVel() {
        return yVel;
    }

}
