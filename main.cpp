#include <iostream>
#include <math.h>
#include <iomanip>
using namespace std;

void getMinor(double **mas, double **p, int i, int j, int n) {
    int ki, kj, di, dj;
    di = 0;
    for (ki = 0; ki < n - 1; ki++) {
        if (ki == i)
            di = 1;
        dj = 0;
        for (kj = 0; kj < n - 1; kj++) {
            if (kj == j)
                dj = 1;
            p[ki][kj] = mas[ki + di][kj + dj];
        }
    }
}

double getDeter(double **mas, int m) {

    int i, j, k, n;
    double d;
    double **p;
    p = new double*[m];
    for (i = 0; i<m; i++)
        p[i] = new double[m];
    j = 0; d = 0;
    k = 1;
    n = m - 1;
    if (m < 1) cout << "Can't calculate determinant " << p[0][0] << endl;
    if (m == 1) {
        return(mas[0][0]);
    }
    if (m == 2) {
        return(mas[0][0] * mas[1][1] - (mas[1][0] * mas[0][1]));
    }

    if (m > 2) {
            for (i = 0; i<m; i++) {
                getMinor(mas, p, i, 0, m);
                d = d + k * mas[i][0] * getDeter(p, n);
                k = -k;
            }
        return(d);
        }
}



double** transp(double **Matrix, int n, int m)
{
    double **Matrix2 = new double*[m];
    for (int i = 0; i < m; i++){
        Matrix2[i] = new double [n];
    }
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            Matrix2[i][j] = Matrix[j][i];

    return Matrix2;
}

double** multiplyMatrx(double **masF, double **masS, int n1, int m1, int n2, int m2){ //m - столбцы в 1,n - количество строк masS, второй матрицы
    double** result = new double*[n1];
    for (int i = 0; i < n1; i++){
        result[i] = new double[n2];
    }

    for (int i = 0; i < n2; i++)
        for (int j = 0; j < m1; j++){
            double c = 0;
            for (int l = 0; l < n2; l++){
                c += masF[i][l] * masS[l][j];
            }
            result[i][j] = c;
    }
    return result;
}

double* multiplyMatrx4d(double **masF, double *masS, int n1, int m1, int n2, int m2){ //m - столбцы в 1,n - количество строк masS, второй матрицы
    double* result = new double[n1];

    for (int i = 0; i < n1; i++){
            double c = 0;
            for (int l = 0; l < n2; l++){
                    c += masF[i][l] * masS[l];
            }
            result[i] = c;
        }
    return result;
}




double** invMatrx(double **matr, int n) {

    long double det;
    det = getDeter(matr, n);

    double **p;
    p = new double *[n];
    double **b;
    b = new double *[n];
    double **c;
    c = new double *[n];
    double **obr_matr = new double *[n];
    double **tobr_matr = new double *[n];
    for (int i = 0; i < n; i++) {
        obr_matr[i] = new double[n];
        tobr_matr[i] = new double[n];
        p[i] = new double[n];
        b[i] = new double[n];
        c[i] = new double[n];
    }
    int k = 1;


    if (det) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                int m = n - 1;

                getMinor(matr, p, i, j, n);

                c[i][j] = k * getDeter(p, m);
                k = -k;
            }
        }

        c = transp(c, n, n);

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                obr_matr[i][j] = c[i][j] / (double) (det);
            }
        }


        obr_matr = transp(obr_matr, n, n);

        return obr_matr;
    }
}

    double *findCoef(double **X, double *d) {
        int n = 3, m = 3;
        double **X2 = new double *[m];
        for (int i = 0; i < m; i++) {
            X2[i] = new double[n];
        }

        double **Xt = new double *[n];
        for (int i = 0; i < n; i++) {
            Xt[i] = new double[m];
        }

        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                Xt[i][j] = X[i][j];

        double **X3 = new double *[m];
        for (int i = 0; i < m; i++) {
            X3[i] = new double[m];
        }

        double **Xi = new double *[m];
        for (int i = 0; i < m; i++) {
            Xi[i] = new double[m];
        }

        double *A = new double[m];

        X2 = transp(X, n, m);

        X3 = multiplyMatrx(X2, X, n, m, m, n);

        X3 = invMatrx(X3, n);

        X3 = multiplyMatrx(X3, X2, n, n, n, m);

        A = multiplyMatrx4d(X3, d, n, m, n, 1);

        return A;
    }

double **findX(double *X, double *Y, int n) {

    double **Xf = new double *[3];
    for (int i = 0; i < 3; i++) {
        Xf[i] = new double[3];
    }

    for (int i = 0; i < 3; i++)
        Xf[i][0] = 1;

    for (int i = 0; i < 3; i++)
            Xf[i][1] = ((X[i + 1]/(X[i + 1] + Y[i + 1])) + (X[i]/(X[i] + Y[i])))/2;

    for (int i = 0; i < 3; i++)
        Xf[i][2] = ((Y[i + 1]/(X[i + 1] + Y[i + 1])) + (Y[i]/(X[i] + Y[i])))/2;

    return Xf;
}


double *findLx(double *X, double *Y, int n) {

    double *L = new double [3];

    for (int i = 0; i < 3; i++)
        L[i] = log((X[i + 1]/(X[i + 1] + Y[i + 1]))) - log((X[i]/(X[i] + Y[i])));

    return L;
}

double *findLy(double *X, double *Y, int n) {

    double *L = new double [3];

    for (int i = 0; i < 3; i++)
        L[i] = log((Y[i + 1]/(X[i + 1] + Y[i + 1]))) - log((Y[i]/(X[i] + Y[i])));

    return L;
}

double fByX(double x, double y, double* A){
    return x * (A[0] + A[1] * x + A[2] * y);
}

double fByY(double x, double y, double* A){
    return y * (A[0] + A[1] * x + A[2] * y);
}

void rgkMethod(){
    double k0h2,k1h2,k2h2,k3h2,k4h2,k5h2,l0h2,l1h2,l2h2,l3h2,l4h2,l5h2;
    double k0,k1,k2,k3,k4,k5,l0,l1,l2,l3,l4,l5;
    double xnh2, ynh2, xn, yn, error, maxErrorByX, maxErrorByY;

    xnh2 = 0.83;
    ynh2 = 0.17;
    xn = 0.83;
    yn = 0.17;
    int n = 3, m = 3;

    double **Xf = new double * [n];
    for(int i = 0; i < n; i++){
        Xf[i] = new double[m];
    }

    double *Xt = new double [n];
    double *Yt = new double [n];
    double *lX = new double [n];
    double *lY = new double [n];

    cout << "Enter data for 1st company: \n";
    for(int i = 0; i < n; i++)
        cin >> Xt[i];

    cout << "Enter data for 2nd company: \n";
    for(int i = 0; i < n; i++)
        cin >> Yt[i];

    /*Xt[0] = 5.35;
    Xt[1] = 5.96;
    Xt[2] = 5.81;
    Xt[3] = 5.18;

    Yt[0] = 0.646;
    Yt[1] = 1.31;
    Yt[2] = 1.75;
    Yt[3] = 1.09;*/
    Xf = findX(Xt, Yt, n);
    lX = findLx(Xt, Yt, n);
    lY = findLy(Xt, Yt, n);



    double *Ax = new double [n];
    double *Ay = new double [n];
    Ax = findCoef(Xf, lX);
    Ay = findCoef(Xf, lY);

    double h = 0.1;
    double h2 = 0.05;

    cout << setprecision(20) << "x0 : " << xnh2 << endl;
    cout << "y0 : " << ynh2 << endl;

    for(int i = 1; i < 11; i++){
        k0h2 = fByX(xnh2, ynh2, Ax);
        l0h2 = fByY(xnh2, ynh2, Ay);

        k1h2 = fByX(xnh2 + h2 * ((double) 1 / 4 * k0h2), ynh2 + h2 * ((double) 1 / 4 * l0h2), Ax);
        l1h2 = fByY(xnh2 + h2 * ((double) 1 / 4 * k0h2), ynh2 + h2 * ((double) 1 / 4 * l0h2), Ay);

        k2h2 = fByX(xnh2 + h2 * ((double) 3 / 40 * k0h2 + (double) 9 / 40 * k1h2),
                  ynh2 + h2 * ((double) 3 / 40 * l0h2 + (double) 9 / 40 * l1h2), Ax);
        l2h2 = fByY(xnh2 + h2 * ((double) 3 / 40 * k0h2 + (double) 9 / 40 * k1h2),
                  ynh2 + h2 * ((double) 3 / 40 * l0h2 + (double) 9 / 40 * l1h2), Ay);

        k3h2 = fByX(xnh2 + h2 * (0.3 * k0h2 - 0.9 * k1h2 + 1.2 * k2h2), ynh2 + h2 * (0.3 * l0h2 - 0.9 * l1h2 + 1.2 * l2h2), Ax);
        l3h2 = fByY(xnh2 + h2 * (0.3 * k0h2 - 0.9 * k1h2 + 1.2 * k2h2), ynh2 + h2 * (0.3 * l0h2 - 0.9 * l1h2 + 1.2 * l2h2), Ay);

        k4h2 = fByX(xnh2 + h2 * ((double) 226 / 729 * k0h2 - (double) 25 / 27 * k1h2 + (double) 880 / 729 * k2h2 +
                               (double) 55 / 729 * k3h2),
                  ynh2 + h2 * ((double) 226 / 729 * l0h2 - (double) 25 / 27 * l1h2 + (double) 880 / 729 * l2h2 +
                               (double) 55 / 729 * l3h2), Ax);
        l4h2 = fByY(xnh2 + h2 * ((double) 226 / 729 * k0h2 - (double) 25 / 27 * k1h2 + (double) 880 / 729 * k2h2 +
                               (double) 55 / 729 * k3h2),
                  ynh2 + h2 * ((double) 226 / 729 * l0h2 - (double) 25 / 27 * l1h2 + (double) 880 / 729 * l2h2 +
                               (double) 55 / 729 * l3h2), Ay);

        k5h2 = fByX(xnh2 + h2 * (-(double) 181 / 270 * k0h2 + (double) 5 / 2 * k1h2 - (double) 266 / 297 * k2h2 -
                               (double) 91 / 27 * k3h2 + (double) 189 / 55 * k4h2),
                  ynh2 + h2 * (-(double) 181 / 270 * k0h2 + (double) 5 / 2 * k1h2 - (double) 266 / 297 * k2h2 -
                               (double) 91 / 27 * k3h2 + (double) 189 / 55 * k4h2), Ax);
        l5h2 = fByY(xnh2 + h2 * (-(double) 181 / 270 * k0h2 + (double) 5 / 2 * k1h2 - (double) 266 / 297 * k2h2 -
                               (double) 91 / 27 * k3h2 + (double) 189 / 55 * k4h2),
                  ynh2 + h2 * (-(double) 181 / 270 * k0h2 + (double) 5 / 2 * k1h2 - (double) 266 / 297 * k2h2 -
                               (double) 91 / 27 * k3h2 + (double) 189 / 55 * k4h2), Ay);

        xnh2 = xnh2 + h2 * ((double)31 / 540 * k0h2 + 0 * k1h2 + (double)190 / 297 * k2h2 - (double)145 / 108 * k3h2 + (double)351 / 220 * k4h2 + (double)1 / 20 * k5h2);
        ynh2 = ynh2 + h2 * ((double)31 / 540 * l0h2 + 0 * l1h2 + (double)190 / 297 * l2h2 - (double)145 / 108 * l3h2 + (double)351 / 220 * l4h2 + (double)1 / 20 * l5h2);

        k0 = fByX(xn, yn, Ax);
        l0 = fByY(xn, yn, Ay);

        k1 = fByX(xn + h * ((double) 1 / 4 * k0), yn + h * ((double) 1 / 4 * l0), Ax);
        l1 = fByY(xn + h * ((double) 1 / 4 * k0), yn + h * ((double) 1 / 4 * l0), Ay);

        k2 = fByX(xn + h * ((double) 3 / 40 * k0 + (double) 9 / 40 * k1),
                  yn + h * ((double) 3 / 40 * l0 + (double) 9 / 40 * l1), Ax);
        l2 = fByY(xn + h * ((double) 3 / 40 * k0 + (double) 9 / 40 * k1),
                  yn + h * ((double) 3 / 40 * l0 + (double) 9 / 40 * l1), Ay);

        k3 = fByX(xn + h * (0.3 * k0 - 0.9 * k1 + 1.2 * k2), yn + h * (0.3 * l0 - 0.9 * l1 + 1.2 * l2), Ax);
        l3 = fByY(xn + h * (0.3 * k0 - 0.9 * k1 + 1.2 * k2), yn + h * (0.3 * l0 - 0.9 * l1 + 1.2 * l2), Ay);

        k4 = fByX(xn + h * ((double) 226 / 729 * k0 - (double) 25 / 27 * k1 + (double) 880 / 729 * k2 +
                            (double) 55 / 729 * k3),
                  yn + h * ((double) 226 / 729 * l0 - (double) 25 / 27 * l1 + (double) 880 / 729 * l2 +
                            (double) 55 / 729 * l3), Ax);
        l4 = fByY(xn + h * ((double) 226 / 729 * k0 - (double) 25 / 27 * k1 + (double) 880 / 729 * k2 +
                            (double) 55 / 729 * k3),
                  yn + h * ((double) 226 / 729 * l0 - (double) 25 / 27 * l1 + (double) 880 / 729 * l2 +
                            (double) 55 / 729 * l3), Ay);

        k5 = fByX(xn + h * (-(double) 181 / 270 * k0 + (double) 5 / 2 * k1 - (double) 266 / 297 * k2 -
                            (double) 91 / 27 * k3 + (double) 189 / 55 * k4),
                  yn + h * (-(double) 181 / 270 * k0 + (double) 5 / 2 * k1 - (double) 266 / 297 * k2 -
                            (double) 91 / 27 * k3 + (double) 189 / 55 * k4), Ax);
        l5 = fByY(xn + h * (-(double) 181 / 270 * k0 + (double) 5 / 2 * k1 - (double) 266 / 297 * k2 -
                            (double) 91 / 27 * k3 + (double) 189 / 55 * k4),
                  yn + h * (-(double) 181 / 270 * k0 + (double) 5 / 2 * k1 - (double) 266 / 297 * k2 -
                            (double) 91 / 27 * k3 + (double) 189 / 55 * k4), Ay);

        xn = xn + h * ((double)31/540 * k0 + 0 * k1 + (double)190/297 * k2 - (double)145/108 * k3 + (double)351/220 * k4 + (double)1/20 * k5);
        yn = yn + h * ((double)31/540 * l0 + 0 * l1 + (double)190/297 * l2 - (double)145/108 * l3 + (double)351/220 * l4 + (double)1/20 * l5);

        cout << "x(" << i << ") : " << xn << endl;
        cout << "y(" << i << ") : " << yn << endl;

        error = abs(xn - xnh2) / (pow(2,9) - 1);
        if (maxErrorByX < error)
            maxErrorByX = error;
        error = abs(yn - ynh2) / (pow(2,9) - 1);
        if (maxErrorByY < error)
            maxErrorByY = error;

    }
    cout << "\nmaxError of x: " << maxErrorByX;
    cout << "\nmaxError of y: " << maxErrorByY;
}
int main() {

    rgkMethod();
    system("pause");
    return 0;
}
