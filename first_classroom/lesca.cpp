//Name: Lesca, Enrico
//Date: 06/11/2026
//Code Output:

//Trapx L (n=10) =7.122589987692e+00
//Trapx L (n=100) =7.535183280883e+00
//Trapx L (n=1000) =7.540331004664e+00
//Gauss L (n=10) =7.344182994898e+00
//Gauss L (n=100) =7.540217122615e+00
//Gauss L (n=1000) =7.540331571809e+00

#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

double eccentricity(double E);
double secant(double (*F)(double),double a,double b,double tol);
double trapezoid(double (*F)(double),double a,double b,int n);
double gauss(double (*F)(double),double a,double b,int n,int ng);
double velocity(double T);

double t; // variabile globale
int main(){
    cout <<setprecision(12);
    cout <<setiosflags(ios:: scientific);
    double t0 =0;
    double thi = 10;
    int n =10;
    int ng = 3;
    double sum1,sum2;
    cout<<"**************************************"<<endl;
    for(int i=10;i<=1000;i*=10){
        sum1 =trapezoid(velocity,t0,thi,i);
        cout <<"Trapx L (n="<<i<<  ") =" <<sum1 <<endl;

    }
    for(int i=10;i<=1000;i*=10){
        sum2=gauss(velocity,t0,thi,i,ng);
        cout <<"Gauss L (n="<<i<<  ") =" <<sum2 <<endl;

    }
    cout <<"*****************************************"<<endl;
    

    return 0;
}

double eccentricity(double E){
    double e = 0.85;
    double alpha = 1;
    return E -e*sin(E) -alpha*t;
}
double secant(double (*F)(double),double a,double b,double tol){
    double fa = F(a);
    double fb = F(b);
    double c = (b+a)/2;
    double fc = F(c);
    double zero;
    //uso bisezione per trovare intervallo migliore per la guess
    if(fabs(fa*fc)<0){
        b =c;
        fb =fc;
    }
    if(fabs(fa*fc)>0){
        a =c;
        fa = fc;
    }
    int n =0;
    double dx = b-a;
    
    while(fabs(dx)>tol && n<100){
        //blocco se dx =0;
        if(fabs(b-a)<1e-14){
            break;
        }
        dx = fb*(b-a)/(fb-fa);
        a =b;
        fa=fb;
        b = b-dx;
        fb =F(b);
        n++;
    }
    zero = b;
    return zero;

}
double velocity(double T){
    t =T;
    double tol = 1e-10;
    double zero =secant(eccentricity,0,t,tol);
    
    return sqrt((1+0.85*cos(zero))/(1-0.85*cos(zero)));
}
double trapezoid(double (*F)(double),double a,double b,int n){
    double dx = fabs(b-a)/(double)n;
    double sum =0;
    for(int i=0;i<n;i++){
        sum+=F(a+i*dx+dx/2)*dx;

    }
    return sum;
}
double gauss(double (*F)(double),double a,double b,int n,int ng){
    double x[ng],w[ng];
    double sum =0.0;
    if(ng == 1){
        x[0] = 0.;
        w[0] = 2.;
    }
    else if(ng == 2){
        x[1] = sqrt(1./3.);
        x[0] = -sqrt(1./3.);
        w[0] = 1.;
        w[1] = 1.;
    }
    else if(ng == 3){
        x[0] = -sqrt(3./5.);
        x[1] = 0.;
        x[2] = sqrt(3./5.);
        w[0] = 5./9.;
        w[1] = 8./9.;
        w[2] = 5./9.;
    }
    else if(ng == 4){
        x[0] = -sqrt( 3./7. - 2./7.*sqrt(6./5.) );
        x[1] = -sqrt( 3./7. + 2./7.*sqrt(6./5.) );
        x[2] = sqrt( 3./7. - 2./7.*sqrt(6./5.) );
        x[3] = sqrt( 3./7. + 2./7.*sqrt(6./5.) );
        w[0] = ( 18. + sqrt(30.) )/( 36. );
        w[1] = ( 18. - sqrt(30.) )/( 36. );
        w[2] = ( 18. + sqrt(30.) )/( 36. );
        w[3] = ( 18. - sqrt(30.) )/( 36. );
    }
    else if(ng == 5){
        x[0] = -1./3.*sqrt( 5. - 2.*sqrt(10./7.) );
        x[1] = -1./3.*sqrt( 5. + 2.*sqrt(10./7.) );
        x[2] = 0.;
        x[3] = 1./3.*sqrt( 5. - 2.*sqrt(10./7.) );
        x[4] = 1./3.*sqrt( 5. + 2.*sqrt(10./7.) );

        w[0] = ( 322. + 13.*sqrt(70.) )/( 900. );
        w[1] = ( 322. - 13.*sqrt(70.) )/( 900. ); 
        w[2] = 128./225.;
        w[3] = ( 322. + 13.*sqrt(70.) )/( 900. );
        w[4] = ( 322. - 13.*sqrt(70.) )/( 900. );
    }
    else{
        cout << "codice non implementato per Ng>5." << endl;
        return 0;
    }
    double dx = fabs(b-a)/(double) n;
    double sumj = 0.0;
    for(int i=0;i<n;i++){
        double x0 = a+i*dx;
        double x1 = x0 +dx;
        double s1 = (x1-x0)/2;
        double s2 = (x1+x0)/2;
        sumj =0.0; //inizializzo la somma a ogni ciclo
        for(int j =0;j<ng;j++){
            sumj += s1*w[j]*F(s1*x[j]+s2);
        }
        sum += sumj;
    }
    return sum;


}

