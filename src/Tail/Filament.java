package Tail;

import ij.process.FloatProcessor;
import ij.process.ImageProcessor;
import java.util.ArrayList;
import java.util.Random;

public class Filament {

    private double x, y, hgu;
    int length;
    ArrayList xPix = new ArrayList();
    ArrayList yPix = new ArrayList();
    ImageProcessor plot;
    Random R = new Random();
    boolean branchX = true;
    private double angle, branch = -70, branchOffset;

    public Filament(double xc, double yc, double a0, ImageProcessor ip, double hgu) {
        this.x = xc;
        this.y = yc;
        this.angle = a0;
        this.plot = ip;
        this.length = 0;
        this.hgu = hgu;
        branchOffset = R.nextGaussian() * hgu * 0.2;
        if (R.nextBoolean()) {
            branchOffset *= -1;
        }
    }

    public void grow(FloatProcessor densityField, FloatProcessor nutrientField, double dfGradSens, double noise, double nfGradSens) {
        xPix.add(new Double(x));
        yPix.add(new Double(y));
        double nfParams[] = getVector(nutrientField, x, y);
//        double dfParams[] = getVector(densityField, x, y);
        if (nfParams == null) {
            nfParams = new double[2];
            nfParams[0] = 0.0;
            nfParams[1] = angle;
        }
//        if (dfParams == null) {
//            dfParams = new double[2];
//            dfParams[0] = 0.0;
//            dfParams[1] = angle;
//        }
//        angle += (dfParams[0] * (dfParams[1] - angle) / dfGradSens)
//                + (nfParams[0] * (nfParams[1] - angle) / nfGradSens)
//                + noise * R.nextGaussian();
        angle += noise * R.nextGaussian();

        double xVec = Math.cos(Math.toRadians(angle));
        double yVec = Math.sin(Math.toRadians(angle));

        x += xVec;
        y += yVec;
        plot.drawLine((int) Math.round(((Double) xPix.get(length)).doubleValue()),
                (int) Math.round(((Double) yPix.get(length)).doubleValue()),
                (int) Math.round(x), (int) Math.round(y));
        length++;
    }

    public int getLength() {
        return length;
    }

    public void resetLength() {
        Double tempX = ((Double) xPix.get(length - 1));
        Double tempY = ((Double) yPix.get(length - 1));
        xPix.clear();
        yPix.clear();
        xPix.add(tempX);
        yPix.add(tempY);
        length = 1;
    }

    public int getBranchx(double d) {
        return (int) Math.round((Double) xPix.get(length - (int) Math.round(d) - 1));
    }

    public int getBranchy(double d) {
        return (int) Math.round((Double) yPix.get(length - (int) Math.round(d) - 1));
    }

    public double getBranchAngle() {
        branch += -2 * branch;
        return (branch + angle);
    }

    public int getX() {
        return (int) x;
    }

    public int getY() {
        return (int) y;
    }

    public double getBranchOffset() {
        return branchOffset;
    }

    public void resetBranchOffset() {
        branchOffset = R.nextGaussian() * hgu * 0.2;
        if (R.nextBoolean()) {
            branchOffset *= -1;
        }
    }

    double[] getVector(FloatProcessor field, double x, double y) {
        double xGrad = field.getInterpolatedValue(x + 1.0, y - 1.0)
                + 2.0 * field.getInterpolatedValue(x + 1.0, y)
                + field.getInterpolatedValue(x + 1.0, y + 1.0)
                - field.getInterpolatedValue(x - 1.0, y - 1.0)
                - 2.0 * field.getInterpolatedValue(x - 1.0, y)
                - field.getInterpolatedValue(x - 1.0, y + 1.0);
        double yGrad = field.getInterpolatedValue(x - 1.0, y + 1.0)
                + 2.0 * field.getInterpolatedValue(x, y + 1.0)
                + field.getInterpolatedValue(x + 1.0, y + 1.0)
                - field.getInterpolatedValue(x - 1.0, y - 1.0)
                - 2.0 * field.getInterpolatedValue(x, y - 1.0)
                - field.getInterpolatedValue(x + 1.0, y - 1.0);
        double gradMag = Math.sqrt(xGrad * xGrad + yGrad * yGrad);

        double t;
        if (gradMag > 0.0) {
            t = Math.toDegrees(Math.atan(yGrad / xGrad));
        } else {
            t = Double.NaN;
        }
        if (xGrad > 0.0 && yGrad < 0.0) {
            t += 360.0;
        } else if (xGrad < 0.0) {
            t += 180.0;
        }
        double params[] = {gradMag, t};
        if (Double.isNaN(t)) {
            params = null;
        }

        return params;
    }
}
