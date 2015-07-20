/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Tail;

import java.util.Random;

/**
 *
 * @author David Barry <david.barry at crick.ac.uk>
 */
public class Virus {

    private double x, y, xVel, yVel;
//    private final double mass = 10e-15; // Mass of single vaccinia virion is ~10 fg - dx.doi.org/10.1016/j.snb.2005.08.047
    private final double mass = 10.0;
    private final double radius = 150.0; // Width of virion is 250 - 350nm
    private Force dragForce;
    private final double rho = 1.0;
    private final double area = Math.PI * radius;
    private final double cD = 0.0005;
    Random r = new Random();

    public Virus(double x, double y) {
        this.x = x;
        this.y = y;
        xVel = 0.0;
        yVel = 0.0;
    }

    public void updateVelocity(Force filamentForce) {
        updateDragForce();
//        System.out.println("xF: " + filamentForce.getxF() + " yF: " + filamentForce.getyF()
//                + " xD: " + dragForce.getxF() + " yD: " + dragForce.getyF()
//                + " xVel: " + xVel + " yVel: " + yVel);
        double xA = r.nextGaussian() + (filamentForce.getxF() + dragForce.getxF()) / mass;
        double yA = r.nextGaussian() + (filamentForce.getyF() + dragForce.getyF()) / mass;
        xVel += xA;
        yVel += yA;
    }

    public void updatePosition() {
        x += xVel;
        y += yVel;
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    public double getRadius() {
        return radius;
    }

    void updateDragForce() {
        double xDF = 0.5 * rho * area * cD * xVel * xVel;
        if (xVel > 0.0) {
            xDF *= -1.0;
        }
        double yDF = 0.5 * rho * area * cD * yVel * yVel;
        if (yVel > 0.0) {
            yDF *= -1.0;
        }
        dragForce = new Force(xDF, yDF);
    }

    public double getxVel() {
        return xVel;
    }

    public void setxVel(double xVel) {
        this.xVel = xVel;
    }

    public double getyVel() {
        return yVel;
    }

    public void setyVel(double yVel) {
        this.yVel = yVel;
    }

}
