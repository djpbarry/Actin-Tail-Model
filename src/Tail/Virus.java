/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package Tail;

/**
 *
 * @author David Barry <david.barry at crick.ac.uk>
 */
public class Virus {

    private double x, y, xVel, yVel;
    private final double mass = 10e-15; // Mass of single vaccinia virion is ~10 fg - dx.doi.org/10.1016/j.snb.2005.08.047
    private final double radius = 150.0; // Width of virion is 250 - 350nm

    public Virus(double x, double y) {
        this.x = x;
        this.y = y;
        xVel = 0.0;
        yVel = 0.0;
    }

    public void updateVelocity(Force force) {
        double xA = force.getxF() / mass;
        double yA = force.getyF() / mass;
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

}
