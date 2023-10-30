import java.util.Scanner;

class MatrixOperations {

    // Function for creating a matrix
    public static double[][] identity(int n) {
        double[][] result = new double[n][n];
        for (int i = 0; i < n; i++) {
            result[i][i] = 1.0;
        }
        return result;
    }

    // Function for subtracting matrix B fom A
    public static double[][] subtract(double[][] A, double[][] B) {
        int rows = A.length;
        int cols = A[0].length;

        if (B.length != rows || B[0].length != cols) {
            throw new IllegalArgumentException("The method is not applicable!”");
        }

        double[][] result = new double[rows][cols];
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i][j] = A[i][j] - B[i][j];
            }
        }
        return result;
    }

    // Function for multiplying matrices A and B
    public static double[][] multiply(double[][] A, double[][] B) {
        int rowsA = A.length;
        int colsA = A[0].length;
        int rowsB = B.length;
        int colsB = B[0].length;

        if (colsA != rowsB) {
            throw new IllegalArgumentException("The method is not applicable!”");
        }

        double[][] result = new double[rowsA][colsB];

        for (int i = 0; i < rowsA; i++) {
            for (int j = 0; j < colsB; j++) {
                result[i][j] = 0;
                for (int k = 0; k < colsA; k++) {
                    result[i][j] += A[i][k] * B[k][j];
                }
            }
        }

        return result;
    }


    // Function for multiplying matrix A and vector b
    public static double[] multiplyMatrixVector(double[][] A, double[] v) {
        int rowsA = A.length;
        int colsA = A[0].length;
        double[] result = new double[rowsA];

        for (int i = 0; i < rowsA; i++) {
            result[i] = 0;
            for (int j = 0; j < colsA; j++) {
                result[i] += A[i][j] * v[j];
            }
        }

        return result;
    }

    // Function for finding inverse of a matrix
    public static double[][] inverse(double[][] A) {
        int n = A.length;
        double[][] augmentedMatrix = new double[n][2 * n];
        double[][] inverse = new double[n][n];

        // Full matrix version [A | I]
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                augmentedMatrix[i][j] = A[i][j];
                augmentedMatrix[i][j + n] = (i == j) ? 1 : 0;
            }
        }

        // Gauss-Jordan method
        for (int i = 0; i < n; i++) {
            // Max element in a column
            int maxRow = i;
            for (int k = i + 1; k < n; k++) {
                if (Math.abs(augmentedMatrix[k][i]) > Math.abs(augmentedMatrix[maxRow][i])) {
                    maxRow = k;
                }
            }

            // Change row to a max element row and visa versa
            double[] temp = augmentedMatrix[i];
            augmentedMatrix[i] = augmentedMatrix[maxRow];
            augmentedMatrix[maxRow] = temp;

            // Normalize row i
            double factor = augmentedMatrix[i][i];
            for (int j = 0; j < 2 * n; j++) {
                augmentedMatrix[i][j] /= factor;
            }

            // Zero columns above diagonal and under it
            for (int k = 0; k < n; k++) {
                if (k != i) {
                    double factor2 = augmentedMatrix[k][i];
                    for (int j = 0; j < 2 * n; j++) {
                        augmentedMatrix[k][j] -= factor2 * augmentedMatrix[i][j];
                    }
                }
            }
        }

        // Get inverse matrix from right hand side of a full one
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                inverse[i][j] = augmentedMatrix[i][j + n];
            }
        }

        return inverse;
    }


    // Function for finding transpose
    public static double[][] transpose(double[][] A) {
        int rows = A.length;
        int cols = A[0].length;
        double[][] result = new double[cols][rows];

        for (int i = 0; i < cols; i++) {
            for (int j = 0; j < rows; j++) {
                result[i][j] = A[j][i];
            }
        }

        return result;
    }

    // Function for subtracting vectors
    public static double[] subtractVectors(double[] a, double[] b) {
        double[] result = new double[a.length];
        for (int i = 0; i < a.length; i++) {
            result[i] = a[i] - b[i];
        }
        return result;
    }

    // Function for creating diagonal matrix
    public static double[][] createDiagonalMatrixFromVector(double[] vector) {
        int n = vector.length;
        double[][] matrix = new double[n][n];
        for (int i = 0; i < n; i++) {
            matrix[i][i] = vector[i];
        }
        return matrix;
    }

}

public class InteriorPointAlgorithm {
    public static double vectorNorm(double[] vector) {
        double sum = 0;
        for (double value : vector) {
            sum += value * value;
        }
        return Math.sqrt(sum);
    }

    public static void main(String[] args) {
        Scanner scanner = new Scanner(System.in);

        // Input for matrix D
        System.out.println("Enter the number of variables in initial trial solution (n):");
        int n = scanner.nextInt();
        double[][] D;
        double[] input = new double[n];
        System.out.println("Enter the initial trial solution:");
        for (int i = 0; i < n; i++) {
            input[i] = scanner.nextDouble();
        }
        D = MatrixOperations.createDiagonalMatrixFromVector(input);

        // Input for matrix A
        System.out.println("Enter the number of inequalities:");
        int m = scanner.nextInt();
        double[][] A = new double[m][n];
        System.out.println("Enter the elements of matrix A:");
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = scanner.nextDouble();
            }
        }

        // Input for vector c
        double[] c = new double[n];
        System.out.println("Enter the elements of vector c:");
        for (int i = 0; i < n; i++) {
            c[i] = scanner.nextDouble();
        }


        solve(A, D, c, 0.5, n);
        System.out.println();
        solve(A, D, c, 0.9, n);

    }

    public static void solve(double[][] A, double[][] D, double[] c, double alpha, int n){
        double[] x = new double[n];
        double[] prevX = new double[n];
        double norm;
        do {
            // Calculate AD
            double[][] AD = MatrixOperations.multiply(A, D);

            // Calculate Dc
            double[] Dc = MatrixOperations.multiplyMatrixVector(D, c);

            // Calculate P
            double[][] ADT = MatrixOperations.transpose(AD);
            double[][] ADADT = MatrixOperations.multiply(AD, ADT);
            double[][] ADADTInv = MatrixOperations.inverse(ADADT);
            double[][] ADTADADTInv = MatrixOperations.multiply(ADT, ADADTInv);
            double[][] ADTADADTInvAD = MatrixOperations.multiply(ADTADADTInv, AD);
            double[][] I = MatrixOperations.identity(n);
            double[][] P = MatrixOperations.subtract(I, ADTADADTInvAD);
            double[]PDc = MatrixOperations.multiplyMatrixVector(P, Dc);

            // Identify the negative component of cp having the largest absolute value
            double v = identifyLargestNegativeComponent(PDc);

            // Compute the x vector
            double alphaV = alpha / v;
            for (int i = 0; i < PDc.length; i++) {
                PDc[i] = PDc[i] * alphaV;
                x[i] = 1 + PDc[i];
            }
            double[] Dx = MatrixOperations.multiplyMatrixVector(D, x);

            double[] diff = MatrixOperations.subtractVectors(Dx, prevX);
            norm = vectorNorm(diff);

            System.arraycopy(Dx, 0, prevX, 0, n);
            D = MatrixOperations.createDiagonalMatrixFromVector(Dx);
        } while (norm >= 0.00001);

        System.out.println("Vector x for alpha " + alpha + ":");
        printVector(prevX);
        calculateAndPrintResult(prevX, c);
    }

    public static void calculateAndPrintResult(double[] vector, double[] c){
        double result = 0;
        for(int i = 0; i < vector.length; i++){
            result += vector[i] * c[i];
        }
        System.out.println("Optimized value of objective function equals " + result);
    }

    public static double identifyLargestNegativeComponent(double[] vector) {
        double largestNegativeValue = 0;
        for (double value : vector) {
            if (largestNegativeValue < Math.abs(value) && value < 0) {
                largestNegativeValue = Math.abs(value);
            }
        }
        return largestNegativeValue;
    }

    public static void printMatrix(double[][] matrix) {
        for (double[] row : matrix) {
            for (double val : row) {
                System.out.print(val + " ");
            }
            System.out.println();
        }
    }

    public static void printVector(double[] vector) {
        for (double val : vector) {
            System.out.println(val);
        }
    }
}
