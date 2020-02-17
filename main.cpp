#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <ctime>
#include <iomanip>
using namespace std;
typedef double (*function) (double);
double DoubleRand (double _max, double _min)
{   return _min + double (rand ()) / RAND_MAX * (_max - _min); }


double fun_X (double y, double z, double a, double b, double c)
{
    return sqrt(1 - pow(y/b, 2)- pow(z/c, 2) ) * a;
}

double fun_Y (double x, double z, double a, double b, double c)
{
    return sqrt(1 - pow(x/a, 2)- pow(z/c, 2) ) * b;
}

double fun_Z (double x, double y, double a, double b, double c)
{
    return sqrt(1 - pow(x/a, 2)- pow(y/b, 2) ) * c;
}

double Probeg(double gamma){
    double Scat = 3.346e-2; //1.599E-05	3.346E-02 coherent and non-coherent respectevly
    double pho = 2.6989;
    return -log(abs(DoubleRand(1,0)))/(Scat*pho);
}

double*** Monte_karlo1 (double a, double b, double c, int N)

{
    auto ***Ellips_Koords = new double** [N]; // N строк в массиве
    for (int count = 0; count < N+1; count++){
        Ellips_Koords[count] = new double* [2];
        for (int count_1 = 0; count_1 < 2; count_1++)
            Ellips_Koords[count][count_1] = new double [3];
    }
    ofstream fout("data_elips.txt");
    fout<<"[";
    double pi = 3.141592653589793;
    int i = 0;
    while (i <= N)
    {      double x = 2 * a * DoubleRand(1, 0) - a;
           double y = 2 * b * DoubleRand(1, 0) - b;
           double z = 2 * c * DoubleRand(1, 0) - c;
           if (((x > 0 and x < fun_X(y, z, a, b, c)) or (x < 0 and x > -fun_X(y, z, a, b, c))) and
               ((z > 0 and z < fun_Z(x, y, a, b, c)) or (z < 0 and z > -fun_Z(x, y, a, b, c))) and
               ((y > 0 and y < fun_Y(x, z, a, b, c)) or (y < 0 and y > -fun_Y(x, z, a, b, c)))) {
               double gamma = DoubleRand(1, 0);
               double gamma_2 = DoubleRand(1, 0);
               cout << gamma << "\t" << gamma_2 << endl;
               if (abs(x + Probeg(gamma) * sin(pi * gamma) * cos(2 * pi * gamma_2)) < b
                   and abs(y + Probeg(gamma) * sin(pi * gamma) * sin(2 * pi * gamma_2)) < a
                   and z + Probeg(gamma) * cos(pi * gamma) > 0) {
                   Ellips_Koords[i][0][0] = x;
                   Ellips_Koords[i][0][1] = y;
                   Ellips_Koords[i][0][2] = z;
                   Ellips_Koords[i][1][0] =
                           Ellips_Koords[i][0][0] + Probeg(gamma) * sin(pi * gamma) * cos(2 * pi * gamma_2);
                   Ellips_Koords[i][1][1] =
                           Ellips_Koords[i][0][1] + Probeg(gamma) * sin(pi * gamma) * sin(2 * pi * gamma_2);
                   Ellips_Koords[i][1][2] = Ellips_Koords[i][0][2] + Probeg(gamma) * cos(pi * gamma);
                   fout << "[";
                   for (int k = 0; k < 2; k++) {
                       fout << "[";
                       for (int j = 0; j < 3; j++) {
                           fout << Ellips_Koords[i][k][j];
                           if (j != 2)
                               fout << ", ";

                       }
                       //if (k!=1)
                       //                        fout<<",";
                       fout << "]";
                       if (k != 1)
                           fout << ", ";
                   }
                   fout << "]";
                   if (i != N)
                       fout << ", ";
                   i = i + 1;
               }
           }
    }
    fout<<"]";
    fout.close(); // закрытие файла
    return Ellips_Koords;
}

int main ()
{   double a = 10; //Parametrs for ellipsoid
    double b = 30;
    double c = 110;
    int N = 200;
    srand(time(0));
    double*** q = Monte_karlo1(a, b, c, N); //cout << abs(4*3.1415*a*b*c/3 - V)/(4*3.1415*a*b*c/3) << endl;
    
}