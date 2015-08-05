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
    private double rE;

    public Energy(double xE, double yE, double rE) {
        this.xE = xE;
        this.yE = yE;
        this.rE = rE;
    }

    public void addEnergy(Energy energy) {
        xE += energy.getxE();
        yE += energy.getyE();
        rE += energy.getwE();
    }

    public double getxE() {
        return xE;
    }

    public double getyE() {
        return yE;
    }

    public double getwE() {
        return rE;
    }

    public double getMag() {
        return Math.sqrt(Math.pow(xE, 2.0) + Math.pow(yE, 2.0));
    }
}
