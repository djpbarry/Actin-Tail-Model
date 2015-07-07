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

    private double x, y, radius, mass, xVel, yVel;

    public Virus(double x, double y, double radius, double mass) {
        this.x = x;
        this.y = y;
        this.radius = radius;
        xVel = 0.0;
        yVel = 0.0;
    }

    public void updateVelocity(double xForce, double yForce) {
        double xA = xForce / mass;
        double yA = yForce / mass;
        xVel += xA;
        yVel += yA;
    }

    public void updatePosition() {
        x += xVel;
        y += yVel;
    }
}
