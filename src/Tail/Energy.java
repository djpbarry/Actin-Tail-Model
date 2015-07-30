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
public class Energy {

    private double xE;
    private double yE;

    public Energy(double xF, double yF) {
        this.xE = xF;
        this.yE = yF;
    }

    public void addEnergy(Energy force) {
        xE += force.getxE();
        yE += force.getyE();
    }

    public double getxE() {
        return xE;
    }

    public double getyE() {
        return yE;
    }
}
