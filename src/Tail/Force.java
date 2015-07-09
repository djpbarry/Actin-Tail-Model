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
public class Force {

    private double xF;
    private double yF;

    public Force(double xF, double yF) {
        this.xF = xF;
        this.yF = yF;
    }

    public void addForce(Force force) {
        xF += force.getxF();
        yF += force.getyF();
    }

    public double getxF() {
        return xF;
    }

    public double getyF() {
        return yF;
    }
}
