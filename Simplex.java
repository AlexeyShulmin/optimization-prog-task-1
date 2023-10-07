import java.util.Scanner;
import java.util.regex.Pattern;
import java.util.regex.Matcher;

public class Simplex {
    private double[][] a;
    private int m, n;

    public Simplex(double[][] A, double[] b, double[] c) {
        m = b.length;
        n = c.length;
        a = new double[m+1][n+m+1];
        for (int i = 0; i < m; i++) {
            for (int j = 0; j < n; j++) {
                a[i][j] = A[i][j];
            }
        }
        for (int i = 0; i < m; i++) a[i][n+i] = 1.0;
        for (int j = 0; j < n; j++) a[m][j] = c[j];
        for (int i = 0; i < m; i++) a[i][m+n] = b[i];
    }

    private int bland() {
        for (int j = 0; j < m+n; j++)
            if (a[m][j] > 0) return j;
        return -1;
    }

    private int minRatioRule(int q) {
        int p = -1;
        for (int i = 0; i < m; i++) {
            if (a[i][q] <= 0) continue;
            else if (p == -1) p = i;
            else if (a[i][m+n] / a[i][q] < a[p][m+n] / a[p][q]) p = i;
        }
        return p;
    }

    public void solve() {
        while (true) {
            int q = bland();
            if (q == -1) break;

            int p = minRatioRule(q);
            if (p == -1) throw new ArithmeticException("Linear program is unbounded");

            pivot(p, q);
        }
    }

    private void pivot(int p, int q) {
        for (int i = 0; i <= m; i++) {
            for (int j = 0; j <= m+n; j++) {
                if (i != p && j != q) a[i][j] -= a[p][j] * a[i][q] / a[p][q];
            }
        }
        for (int i = 0; i <= m; i++)
            if (i != p) a[i][q] = 0.0;
        for (int j = 0; j <= m+n; j++)
            if (j != q) a[p][j] /= a[p][q];
        a[p][q] = 1.0;
    }

    public double value() {
        return -a[m][m+n];
    }

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);

        System.out.println("Enter all the variables in a row in such format: x1, x2, x3, x4");
        String variablesInput = sc.nextLine();
        String[] variables = variablesInput.split(", ");

        System.out.println("Enter your objective function in the following format: z = 3x1 + 4x2 - x3");
        String objectiveFunction = sc.nextLine();

        System.out.println("Is the problem a maximization or minimization? Enter 'max' or 'min':");
        String optimizationType = sc.nextLine();

        System.out.println("Enter the number of inequalities:");
        int numOfInequalities = sc.nextInt();
        sc.nextLine(); // consume the newline

        String[] inequalities = new String[numOfInequalities];
        System.out.println("Enter the inequalities in the following format, one by one: -x1 + 3x2 <= 4");
        for (int i = 0; i < numOfInequalities; i++) {
            inequalities[i] = sc.nextLine();
        }

        // Extract data from user inputs
        double[] c = parseCoefficients(objectiveFunction, variables);
        double[][] A = new double[numOfInequalities][variables.length];
        double[] b = new double[numOfInequalities];

        for (int i = 0; i < numOfInequalities; i++) {
            A[i] = parseCoefficients(inequalities[i], variables);
            b[i] = parseBValue(inequalities[i]);
        }

        if(optimizationType.equals("min")){
            for(int i = 0; i < c.length; i ++){
                c[i] *= -1;
            }
        }

        // Solve the linear program
        Simplex simplex = new Simplex(A, b, c);
        simplex.solve();

        double result = simplex.value();
        System.out.println("Value of objective function: " + String.format("%.2f", result));
    }


    public static double[] parseCoefficients(String s, String[] variables) {
        double[] coefficients = new double[variables.length];
        for (int i = 0; i < variables.length; i++) {
            String regex = "([-+]?\\s?\\d*\\.?\\d*)\\s*\\Q" + variables[i] + "\\E";
            Matcher matcher = Pattern.compile(regex).matcher(s);
            if (matcher.find()) {
                String coefStr = matcher.group(1).replace(" ", "");
                if (coefStr.isEmpty() || "+".equals(coefStr)) {
                    coefficients[i] = 1.0;
                } else if ("-".equals(coefStr)) {
                    coefficients[i] = -1.0;
                } else {
                    coefficients[i] = Double.parseDouble(coefStr);
                }
                if (s.contains(">")) coefficients[i] *= -1;
            } else {
                coefficients[i] = 0.0;
            }
        }
        return coefficients;
    }



    public static double parseBValue(String s) {
        String regex = "=\\s*(-?\\d+)";
        Matcher matcher = Pattern.compile(regex).matcher(s);
        if (matcher.find()) {
            if (s.contains(">"))
                return -1 * Double.parseDouble(matcher.group(1));
            return Double.parseDouble(matcher.group(1));
        }
        return 0; 
    }

    public static void printVector(String name, double[] vec) {
        System.out.print(name + " ");
        for (double v : vec) {
            System.out.print(v + " ");
        }
        System.out.println();
    }

    public static void printMatrix(String name, double[][] mat) {
        System.out.println(name);
        for (double[] row : mat) {
            for (double v : row) {
                System.out.print(v + " ");
            }
            System.out.println();
        }
    }
}
