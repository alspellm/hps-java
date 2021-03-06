package org.hps.recon.tracking.kalman;

/**
 * Simple matrix package strictly for N by N matrices needed by the Kalman fitter. Most of its usage has
 * been replaced by non-object-oriented EJML routines, for speed.
 */
class SquareMatrix { // 
    double[][] M = null;
    int N;

    SquareMatrix(int N) {
        M = new double[N][N];
        this.N = N;
    }

    // Careful here: a new array is not made
    SquareMatrix(int N, double[][] m) {
        M = m;
        this.N = N;
    }

    SquareMatrix(int N, double v) { // Create a diagonal matrix proportional to the unit matrix
        this.N = N;
        M = new double[N][N];
        for (int i = 0; i < N; i++) {
            M[i][i] = v;
        }
    }

    void scale(double f) { // Multiply all matrix elements by a scalar
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                M[i][j] *= f;
            }
        }
    }

    boolean isNaN() {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (Double.isNaN(M[i][j]))
                    return true;
            }
        }
        return false;
    }
    
    SquareMatrix multiply(SquareMatrix M2) { // Standard matrix multiplication
        SquareMatrix Mp = new SquareMatrix(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int k = 0; k < N; k++) {
                    Mp.M[i][j] += M[i][k] * M2.M[k][j];
                }
            }
        }
        return Mp;
    }

    SquareMatrix sum(SquareMatrix m2) { // Add two matrices
        SquareMatrix Ms = new SquareMatrix(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                Ms.M[i][j] = m2.M[i][j] + M[i][j];
            }
        }
        return Ms;
    }

    SquareMatrix dif(SquareMatrix m2) { // Subtract two matrices
        SquareMatrix Ms = new SquareMatrix(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                Ms.M[i][j] = M[i][j] - m2.M[i][j];
            }
        }
        return Ms;
    }

    SquareMatrix transpose() {
        SquareMatrix Mt = new SquareMatrix(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                Mt.M[i][j] = M[j][i];
            }
        }
        return Mt;
    }

    SquareMatrix similarity(SquareMatrix F) { // F * M * F-transpose
        SquareMatrix Mp = new SquareMatrix(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int m = 0; m < N; m++) {
                    for (int n = 0; n < N; n++) {
                        Mp.M[i][j] += F.M[i][m] * M[m][n] * F.M[j][n];
                    }
                }
            }
        }
        return Mp;
    }

    SquareMatrix rotate(RotMatrix R) { // Similarity transform by rotation matrix F
        if (N != 3) {
            System.out.format("SquareMatrix.rotate: incorrect dimension %d\n", N);
            return null;
        }
        SquareMatrix Mp = new SquareMatrix(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int m = 0; m < N; m++) {
                    for (int n = 0; n < N; n++) {
                        Mp.M[i][j] += R.M[i][m] * M[m][n] * R.M[j][n];
                    }
                }
            }
        }
        return Mp;
    }

    SquareMatrix inverseRotate(RotMatrix R) { // Similarity transform by inverted rotation matrix F
        if (N != 3) {
            System.out.format("SquareMatrix.inverseRotate: incorrect dimension %d\n", N);
            return null;
        }
        SquareMatrix Mp = new SquareMatrix(N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                for (int m = 0; m < N; m++) {
                    for (int n = 0; n < N; n++) {
                        Mp.M[i][j] += R.M[m][i] * M[m][n] * R.M[n][j];
                    }
                }
            }
        }
        return Mp;
    }

    void print(String s) {
        System.out.format("%s",this.toString(s));
    }
    
    String toString(String s) {
        String str;
        str = String.format("Printout of matrix %s  %d\n", s, N);
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                str=str+String.format("  %12.4e", M[i][j]);
            }
            str=str+String.format("\n");
        }
        return str;
    }

    SquareMatrix copy() {
        double [][] newArray = new double [N][];
        for (int i=0; i<N; ++i) {
            newArray[i] = new double[N];
            System.arraycopy(M[i], 0, newArray[i], 0, N);  // Doesn't really give any advantage for 5 by 5 arrays
        }
        return new SquareMatrix(N,newArray);
//        SquareMatrix Mc = new SquareMatrix(N);
//        for (int i = 0; i < N; i++) {
//            for (int j = 0; j < N; j++) {
//                Mc.M[i][j] = M[i][j];
//            }
//        }
//        Mc.N = N;
//        return Mc;
    }

    SquareMatrix invert() { // Inversion algorithm copied from "Numerical Methods in C"
        SquareMatrix Minv = this.copy(); // Creates a new matrix; does not overwrite the original one
        return Minv.fastInvert();
    }

    SquareMatrix fastInvert( ) {
        int icol = 0, irow = 0;
        int[] indxc = new int[N];
        int[] indxr = new int[N];
        int[] ipiv = new int[N];
        for (int i = 0; i < N; i++) {
            double big = 0.0;
            for (int j = 0; j < N; j++) {
                if (ipiv[j] != 1) {
                    for (int k = 0; k < N; k++) {
                        if (ipiv[k] == 0) {
                            if (Math.abs(this.M[j][k]) >= big) {
                                big = Math.abs(this.M[j][k]);
                                irow = j;
                                icol = k;
                            }
                        }
                    }
                }

            }
            ++(ipiv[icol]);
            if (irow != icol) {
                for (int l = 0; l < N; l++) {
                    double c = this.M[irow][l];
                    this.M[irow][l] = this.M[icol][l];
                    this.M[icol][l] = c;
                }
            }
            indxr[i] = irow;
            indxc[i] = icol;
            if (this.M[icol][icol] == 0.0) {
                //System.out.format("Singular matrix in SquareMatrix.java method invert.\n");
                return null;
            }
            double pivinv = 1.0 / this.M[icol][icol];
            this.M[icol][icol] = 1.0;
            for (int l = 0; l < N; l++) {
                this.M[icol][l] *= pivinv;
            }
            for (int ll = 0; ll < N; ll++) {
                if (ll != icol) {
                    double dum = this.M[ll][icol];
                    this.M[ll][icol] = 0.0;
                    for (int l = 0; l < N; l++) {
                        this.M[ll][l] -= this.M[icol][l] * dum;
                    }
                }
            }
        }
        for (int l = N - 1; l >= 0; l--) {
            if (indxr[l] != indxc[l]) {
                for (int k = 0; k < N; k++) {
                    double dum = this.M[k][indxr[l]];
                    this.M[k][indxr[l]] = this.M[k][indxc[l]];
                    this.M[k][indxc[l]] = dum;
                }
            }
        }

        return this;
    }
}
