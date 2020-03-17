#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <vector>
#include <ctime>
#include <random>
using namespace std;
std::random_device rd;
std::mt19937 gen(rd());
std::uniform_real_distribution<> urd(0, 1);
<<<<<<< HEAD
double DoubleRand (double _max, double _min){//distribution function
    return _min + urd(gen)*(_max-_min);
}

double fun_ellips (double y, double b, double z, double c, double a){//is point in ellipsoid?
    return sqrt(1 - pow(y/b, 2)- pow(z/c, 2) ) * a;
}

double fun_cirk(double x, double b){//is point in cillinder
    return sqrt(b*b-x*x);
}

double Scat_diff(double alpha, double alpha_sh){//func fo energy
=======
double DoubleRand (double _max, double _min){
    return _min + urd(gen)*(_max-_min);
}

double fun_ellips (double y, double b, double z, double c, double a){
    return sqrt(1 - pow(y/b, 2)- pow(z/c, 2) ) * a;
}

double fun_cirk(double x, double b){
    return sqrt(b*b-x*x);
}

double Scat_diff(double alpha, double alpha_sh){
>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33
    double phi_rcos = 1 - 1 / alpha + 1 / alpha_sh;
    double Sech_d = 0.5 * pow(1 + alpha_sh * (1 - phi_rcos), -2) * (1 + pow(phi_rcos, 2) + (pow(alpha_sh*(1 - phi_rcos), 2))
            / (1 + alpha_sh * (1 - phi_rcos)));
    return 13*0.511*0.511/2*Sech_d;
}

double Scat_int(double alpha){
    double pi = 3.141592653589793;
    return 2*pi*13*0.511*0.511* (1+alpha/(alpha*alpha)*(2*(1+alpha)/(1+2*alpha)-log(1+2*alpha)/alpha)+log(1+2*alpha)
    /(2*alpha)-(1+3*alpha)/pow(1+2*alpha, 2)) ;
<<<<<<< HEAD
}

double kompron_scattering(double energy){
    return -0.01061544616*pow(energy, 3)+0.0658795422219*pow(energy, 2)-0.1388152190602*pow(energy, 1)+0.1431159775481;
}

double total_scat(double energy){
    return 0.050587765 * pow(energy, -0.901791187);
}

double rasst_d(double a,double b,double c){//func return distance between point and detector
=======
}

double kompron_scattering(double energy){
    return -0.01061544616*pow(energy, 3)+0.0658795422219*pow(energy, 2)-0.1388152190602*pow(energy, 1)+0.1431159775481;
}

double total_scat(double energy){
    return 0.050587765 * pow(energy, -0.901791187);
}

double rasst_d(double a,double b,double c){
>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33
    double rast = pow(a,2)+pow(b,2)+pow(c,2);
    while(rast < 1){
        a-=0.5, b-=0.5, c-=0.5;
        rast = pow(a,2)+pow(b,2)+pow(c,2);
    }
    return rast;
}


double Probeg(double gamma, double energy){
    double pho = 2.6989;
    return -log(abs(gamma))/(total_scat(energy)*pho);
}

double Monte_karlo1 (double a, double b, double c, double high_cill, int N, double Start_energy){
    vector<vector<vector<double>>> Ellips_Koords;
    Ellips_Koords.resize(N);
    for (auto & Ellips_Koord : Ellips_Koords) {
        Ellips_Koord.resize(100);
        for (auto & j : Ellips_Koord){
            j.resize(3);
        }
    }
    ofstream fout1("data_elips1.txt");

    double pi = 3.141592653589793;
    int i = 0;
<<<<<<< HEAD
    while (i < N-1){//number of points
=======
    while (i < N-1){
>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33
        vector<vector<double>> Scat_dist;
        Scat_dist.resize(30);
        for (auto & Scat_dist : Scat_dist) {
            Scat_dist.resize(20);
        }
        double x = 2 * a * DoubleRand(1, 0) - a, y = 2 * b * DoubleRand(1, 0) - b;
        double z = 2 * c * DoubleRand(1, 0) - c;
        if (((x > 0 and x < fun_ellips(y, b, z, c, a)) or (x < 0 and x > -fun_ellips(y, b, z, c, a))) and
            ((z > 0 and z < fun_ellips(x, a, y, b, c))) and
<<<<<<< HEAD
            ((y > 0 and y < fun_ellips(x, a, z, c, b)) or (y < 0 and y > -fun_ellips(x, a, z, c, b)))) {//in ellipse?
=======
            ((y > 0 and y < fun_ellips(x, a, z, c, b)) or (y < 0 and y > -fun_ellips(x, a, z, c, b)))) {
>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33
            double alpha = Start_energy/0.511;
            Ellips_Koords[i][0][0] = x;
            Ellips_Koords[i][0][1] = y;
            Ellips_Koords[i][0][2] = z;
<<<<<<< HEAD
            double gamma = DoubleRand(1, 0), gamma_1_probeg = DoubleRand(1, 0);//random values for angles and probeg
=======
            double gamma = DoubleRand(1, 0), gamma_1_probeg = DoubleRand(1, 0);
>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33
            double gamma_2 = DoubleRand(1, 0);
            double x_new = x + Probeg(gamma_1_probeg, Start_energy) * sin(pi * gamma) * cos(2 * pi * gamma_2);
            double y_new = y + Probeg(gamma_1_probeg, Start_energy) * sin(pi * gamma) * sin(2 * pi * gamma_2);
            double z_new = z + Probeg(gamma_1_probeg, Start_energy) * cos(pi * gamma);
<<<<<<< HEAD
            if (abs(x_new) < fun_cirk(y_new, b) and abs(y_new) < fun_cirk(x_new, b) and z_new > 0 and z_new < high_cill) {//in cillinder?
                Ellips_Koords[i][1][0] = x_new;
                Ellips_Koords[i][1][1] = y_new;
                Ellips_Koords[i][1][2] = z_new;
				double W = 1;//start weight
				int j = 0;
                while (alpha > 0.1/0.511) {//until energy of point more then 0.01 mEv
                    if (abs(Ellips_Koords[i][j + 1][0]) < fun_cirk(Ellips_Koords[i][j + 1][1], b)
                    and abs(Ellips_Koords[i][j + 1][1]) < fun_cirk(Ellips_Koords[i][j + 1][0], b)
                        and Ellips_Koords[i][j + 1][2] > 0 and Ellips_Koords[i][j + 1][2] < high_cill) {//in cillinder?
                        double gamma_2_probeg = DoubleRand(1, 0);//random values for angles and probeg
                        double gamma_3 = DoubleRand(1, 0), gamma_4 = DoubleRand(1, 0);
                        if (kompron_scattering(alpha * 0.511) / total_scat(alpha * 0.511) > gamma_3) {//is it kompron scattering?
=======
            if (abs(x_new) < fun_cirk(y_new, b) and abs(y_new) < fun_cirk(x_new, b) and z_new > 0 and z_new < high_cill) {
                Ellips_Koords[i][1][0] = x_new;
                Ellips_Koords[i][1][1] = y_new;
                Ellips_Koords[i][1][2] = z_new;
				double W = 1;
				int j = 0;
                while (alpha > 10/0.511) {
                    if (abs(Ellips_Koords[i][j + 1][0]) < fun_cirk(Ellips_Koords[i][j + 1][1], b)
                    and abs(Ellips_Koords[i][j + 1][1]) < fun_cirk(Ellips_Koords[i][j + 1][0], b)
                        and Ellips_Koords[i][j + 1][2] > 0 and Ellips_Koords[i][j + 1][2] < high_cill) {
                        double gamma_2_probeg = DoubleRand(1, 0);
                        double gamma_3 = DoubleRand(1, 0), gamma_4 = DoubleRand(1, 0);
                        if (kompron_scattering(alpha * 0.511) / total_scat(alpha * 0.511) > gamma_3) {
>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33
                            double alpha_sh = alpha;
                            double rand_1 = DoubleRand(1, 0), rand_2 = DoubleRand(1, 0);
                            alpha = alpha_sh * (1 + 2 * alpha_sh * rand_1) / (1 + 2 * alpha_sh);
                            double p = alpha / alpha_sh + alpha_sh / alpha +
                                       (1 / alpha_sh - 1 / alpha) * (2 + 1 / alpha_sh - 1 / alpha);
<<<<<<< HEAD
                            while (rand_2 * (1 + 2 * alpha_sh + 1 / (1 + 2 * alpha_sh)) >= p) {//energy from description
=======
                            while (rand_2 * (1 + 2 * alpha_sh + 1 / (1 + 2 * alpha_sh)) >= p) {
>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33
                                rand_1 = DoubleRand(1, 0), rand_2 = DoubleRand(1, 0);
                                alpha = alpha_sh * (1 + 2 * alpha_sh * rand_1) / (1 + 2 * alpha_sh);
                                p = alpha / alpha_sh + alpha_sh / alpha +
                                    (1 / alpha_sh - 1 / alpha) * (2 + 1 / alpha_sh - 1 / alpha);
                            }
                            double minus_or_nor = pow(-1, rand());
                            Ellips_Koords[i][j+2][0] =
                                    Ellips_Koords[i][j+1][0] + Probeg(gamma_2_probeg, alpha) * minus_or_nor *
                                                             sin(acos(1 - 1 / alpha + 1 / alpha_sh)) *
                                                             cos(2 * pi * gamma_4);
                            Ellips_Koords[i][j + 2][1] =
                                    Ellips_Koords[i][j + 1][1] + Probeg(gamma_2_probeg, alpha) * minus_or_nor *
                                                                 sin(acos(1 - 1 / alpha + 1 / alpha_sh)) *
                                                                 sin(2 * pi * gamma_4);
                            Ellips_Koords[i][j + 2][2] = Ellips_Koords[i][j + 1][2] +
                                                         Probeg(gamma_2_probeg, alpha) * (1 - 1 / alpha + 1 / alpha_sh);
                            if (abs(Ellips_Koords[i][j + 2][0]) < fun_cirk(Ellips_Koords[i][j + 2][1], b)
                                and abs(Ellips_Koords[i][j + 2][1]) < fun_cirk(Ellips_Koords[i][j + 2][0], b)
<<<<<<< HEAD
                                and Ellips_Koords[i][j + 2][2] > 0 and Ellips_Koords[i][j + 2][2] < high_cill) {//in cillinder, again?
                                W = W*Scat_diff(alpha, alpha_sh)/Scat_int(alpha);//new weight
                                Scat_dist[j+1][0]=alpha*0.511;//new energy
                                Scat_dist[j+1][1]=W*gamma_2_probeg/(rasst_d(Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],Ellips_Koords[i][j + 2][2]));//\nu for differente detectors
										//rasst_d return distance between point and detector
=======
                                and Ellips_Koords[i][j + 2][2] > 0 and Ellips_Koords[i][j + 2][2] < high_cill) {
                                W = W*Scat_diff(alpha, alpha_sh)/Scat_int(alpha);
                                Scat_dist[j+1][0]=alpha*0.511;
                                Scat_dist[j+1][1]=W*gamma_2_probeg/(rasst_d(Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],Ellips_Koords[i][j + 2][2]));
>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33
                                Scat_dist[j+1][2]=W*gamma_2_probeg/(rasst_d(Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill/4 - Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][3]=W*gamma_2_probeg/(rasst_d(Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill/2 - Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][4]=W*gamma_2_probeg/(rasst_d(Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill*3/4 - Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][5]=W*gamma_2_probeg/(rasst_d(Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill-Ellips_Koords[i][j + 2][2]));
<<<<<<< HEAD

                                Scat_dist[j+1][6]=W*gamma_2_probeg/(rasst_d(b-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][7]=W*gamma_2_probeg/(rasst_d(b-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill/4 -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][8]=W*gamma_2_probeg/(rasst_d(b-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill/2 -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][9]=W*gamma_2_probeg/(rasst_d(b-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill*3/4 -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][10]=W*gamma_2_probeg/(rasst_d(b-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill -Ellips_Koords[i][j + 2][2]));

                                Scat_dist[j+1][11]=W*gamma_2_probeg/(rasst_d(b-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][12]=W*gamma_2_probeg/(rasst_d(b/4-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][13]=W*gamma_2_probeg/(rasst_d(b/2 - Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][14]=W*gamma_2_probeg/(rasst_d(b*3/4-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][15]=W*gamma_2_probeg/(rasst_d(Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill -Ellips_Koords[i][j + 2][2]));
                                for(int y=0;y<15;++y){//output in file
                                    if (Scat_dist[j+1][y]!=0.0) {//round energy values
                                        if(y==0)
                                            fout1 <<round( Scat_dist[j + 1][y]*100)/100 << "\t";
                                        else
                                            fout1 << Scat_dist[j + 1][y] << "\t";
                                    }
                                }
                                fout1<<endl;
                                j++;
                                }
                        } else {
                            i = i-1;
                            break;
                        }

=======

                                Scat_dist[j+1][6]=W*gamma_2_probeg/(rasst_d(b-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][7]=W*gamma_2_probeg/(rasst_d(b-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill/4 -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][8]=W*gamma_2_probeg/(rasst_d(b-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill/2 -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][9]=W*gamma_2_probeg/(rasst_d(b-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill*3/4 -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][10]=W*gamma_2_probeg/(rasst_d(b-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill -Ellips_Koords[i][j + 2][2]));

                                Scat_dist[j+1][11]=W*gamma_2_probeg/(rasst_d(b-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][12]=W*gamma_2_probeg/(rasst_d(b/4-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][13]=W*gamma_2_probeg/(rasst_d(b/2 - Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][14]=W*gamma_2_probeg/(rasst_d(b*3/4-Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill -Ellips_Koords[i][j + 2][2]));
                                Scat_dist[j+1][15]=W*gamma_2_probeg/(rasst_d(Ellips_Koords[i][j + 2][0],
                                        Ellips_Koords[i][j + 2][1],high_cill -Ellips_Koords[i][j + 2][2]));
                                for(int y=0;y<15;++y){
                                    if (Scat_dist[j+1][y]!=0.0) {
                                        if(y==0)
                                            fout1 <<round( Scat_dist[j + 1][y]*100)/100 << "\t";
                                        else
                                            fout1 << Scat_dist[j + 1][y] << "\t";
                                    }
                                }
                                fout1<<endl;
                                j++;
                                }
                        } else {
                            i = i-1;
                            break;
                        }

>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33
                    } else {
                        i = i-1;
                        break;
                    }
                }
                i = i + 1;
            }
        }
<<<<<<< HEAD
    }
    ofstream fout("data_elips.txt");
    for(int y=0;y<Ellips_Koords.size()-1;++y){
        for(int x=0;x<Ellips_Koords[y].size();++x){
            for(int z=0;z<Ellips_Koords[y][x].size();++z){
                if (Ellips_Koords[y][x][z]!=0)
                    fout<<Ellips_Koords[y][x][z]<<" ";
            }
        }
        if (y<=Ellips_Koords.size()-1)
        fout<<"|";
    }
    fout.close(); // close
=======
    }
    ofstream fout("data_elips.txt");
    for(int y=0;y<Ellips_Koords.size()-1;++y){
        for(int x=0;x<Ellips_Koords[y].size();++x){
            for(int z=0;z<Ellips_Koords[y][x].size();++z){
                if (Ellips_Koords[y][x][z]!=0)
                    fout<<Ellips_Koords[y][x][z]<<" ";
            }
        }
        if (y<=Ellips_Koords.size()-1)
        fout<<"|";
    }
    fout.close(); // закрытие файла
>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33
    return 0;
}

int main (){
    double a = 10; //Parametrs for ellipsoid
    double b = 30;
    double c = 1e-6;
<<<<<<< HEAD
    double high = 20; //high of cillinder
    double Start_energy = 3.0; //Start energy
    int N = 100000;//number of point in ellipse
    srand(time(0));
    Monte_karlo1(a, b, c, high, N, Start_energy);
}
=======
    double high = 20;
    double Start_energy = 3.0;
    int N = 500;
    srand(time(0));
    Monte_karlo1(a, b, c, high, N, Start_energy);
}
>>>>>>> d082c9b27f01d1aed7a45394379fe92a1fd38c33
