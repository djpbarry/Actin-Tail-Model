/*
 * Copyright (C) 2016 David Barry <david.barry at crick.ac.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package Tail;

import java.util.Random;

/**
 *
 * @author David Barry <david.barry at crick.ac.uk>
 */
public class Monomer {

    public static final int ADP = 0, ATPPI = 1, ATP = 2;
    private final double x, y;
    public static final double DIAMETER = 7;
    private int state;
    private int startTime;
    private final int stateChangeCoeff = 50000;
    private final Random r;

    public Monomer(double x, double y, int startTime) {
        this.x = x;
        this.y = y;
        this.state = Monomer.ADP;
        this.startTime = startTime;
        this.r = new Random();
    }

    public double getX() {
        return x;
    }

    public double getY() {
        return y;
    }

    void changeState(int time) {
        int age = time - startTime;
        if (age > r.nextInt(stateChangeCoeff)) {
            state++;
        }
        this.startTime = 0;
    }

    public int getState() {
        return state;
    }

}
