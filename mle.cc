/*
This program is meant to compute maximal Lyapunov expoennt of Lorenz system.
*/
#include <iostream>
#include <cmath>
#include <complex>
#include <fstream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <blitz/array.h>
#include <sys/stat.h>

using namespace blitz;
// ============ Declaration of functions ===================
double Euclidean_measure(Array<double,1>);
Array<double,1> Lorenz_system(Array<double,1>);
Array<double,1> RK4_integrator(Array<double,1>, double, Array<double,1> (*F)(Array<double, 1>) );
Array<double,1> linearized_Lorenz_system(Array<double,1>, Array<double,1>);
Array<double,1> error_RK4_integrator(Array<double,1>, Array<double,1>, double, Array<double,1> (*F)(Array<double,1>, Array<double,1>));
double max_lypo_expo(Array<double,1> , double, Array<double,1> (*F)(Array<double,1>) );

// =========== I/O files ==================================
ofstream time_series_file;


int main(int argc, char const *argv[])
{
	Array<double,1> X(3);								// State vector for the reference trajectory
	X = -3.16,-5.31,13.31;							// Initial condition r = 28
	// X = 0.62225717,-0.08232857,30.60845379;				// Initial condition r = 45.92
	// X = 4.29033, -0.577049, 28.8962;
	double dt = 0.001;
	max_lypo_expo(X, dt, Lorenz_system);
	return 0;
}

double Euclidean_measure(Array<double,1> A)
{
	int N = A.numElements();
	double d = 0;
	for (int i = 0; i < N; ++i)
	{
		d += pow(A(i),2);
	}
	d = sqrt(d);
	return d;
}

Array<double,1> RK4_integrator(Array<double,1> X0, double dt, Array<double,1> (*F)(Array<double,1>) )
{
	int dim = X0.numElements();
	Array<double,1> X(dim);				// State vector
	// -------- Intermediate derivatives --------
	Array<double,1> K1(dim);
	Array<double,1> K2(dim);
	Array<double,1> K3(dim);
	Array<double,1> K4(dim);
	// ------------------------------------------
	K1 = dt*F(X0);
	X = X0 + 0.5*K1;
	K2 = dt*F(X);
	X = X0 + 0.5*K2;
	K3 = dt*F(X);
	X = X0 + K3;
	K4 = dt*F(X);
	X = X0 + (1.0/6)*(K1 + 2*K2 + 2*K3 + K4);
	return X;

}

Array<double,1> error_RK4_integrator( Array<double,1> E0, Array<double,1> X, double dt, Array<double,1> (*F)(Array<double,1>, Array<double,1>))
{
	// X is the state vector
	int dim = E0.numElements();
	Array<double,1> E(dim);				// Tangent vector
	// -------- Intermediate derivatives --------
	Array<double,1> K1(dim);
	Array<double,1> K2(dim);
	Array<double,1> K3(dim);
	Array<double,1> K4(dim);
	// ------------------------------------------
	K1 = dt*F(E0, X);
	E = E0 + 0.5*K1;
	K2 = dt*F(E0, X);
	E = E0 + 0.5*K2;
	K3 = dt*F(E0, X);
	E = E0 + K3;
	K4 = dt*F(E0, X);
	E = E0 + (1.0/6)*(K1 + 2*K2 + 2*K3 + K4);
	return E;

}

Array<double,1> Lorenz_system(Array<double,1> X)
{
	// ------ Parameters ------------------------
	double sigma = 10;//16; //10;
	double r = 28;//45.92; //28;
	double b = 8.0/3;//4; //8.0/3;
	// ------------------------------------------

	Array<double,1> F(3);
	F(0) = sigma*(-X(0) + X(1));
	F(1) = (r- X(2))*X(0) - X(1);
	F(2) = X(0)*X(1) - b*X(2);
	return F;
}

/*Array<double,1> Rotating_Lorenz_system(Array<double,1> X)
{
	// ------ Parameters ------------------------
	double sigma = 16; //10;
	double r = 45.92; //28;
	double b = 4; //8.0/3;
	// ------------------------------------------

	Array<double,1> F(3);
	F(0) = sigma*(-X(0) + X(1));
	F(1) = (r- X(2))*X(0) - X(1);
	F(2) = X(0)*X(1) - b*X(2);
	return F;
}*/

Array<double,1> linearized_Lorenz_system(Array<double,1> E, Array<double,1> X)
{
	// X is the reference trajectory while
	// E is the error along the refrence trajectory.
	// ------ Parameters ------------------------
	double sigma = 10;
	double r = 28;
	double b = 8.0/3;
	// ------------------------------------------

	Array<double,1> F(3);
	F(0) = sigma*(-E(0) + E(1));
	F(1) = (r- X(2))*E(0) - E(1)- X(0)*E(2);
	F(2) = X(1)*E(0) + X(0)*E(1)- b*E(2);
	return F;
}

double max_lypo_expo(Array<double,1> X, double dt, Array<double,1> (*F)(Array<double,1>) )
{
	
	// Array<double,1> X(3);				// State vector for the reference trajectory
	// X = -3.16,-5.31,13.31;				// Initial condition r = 28
	// X = 0.62225717,-0.08232857,30.60845379;				// Initial condition r = 45.92
	// X = 4.29033, -0.577049, 28.8962;


	int T_transient = 32000;
	// ======  Removing the transience to enter the attractor========
	for (int i = 0; i < T_transient; ++i)
	{
		X = RK4_integrator(X, dt, F);
	}
	// ==============================================================

	cout<<X<<endl;						// After transient, the initial condition to start with.
	// ======================== Calculation of Lyapunov exponent ===============================
	Array<double,1> Y(3);				// State vector for the second trajectory

	Array<double,1> E(3);				// Tangent vector
	E = 1e-8,0,0;						// Initial tangent vector

	double d0;
	d0 = Euclidean_measure(E); 			// Initial error					
	double d = 0;						// length of tangent vector after one iteration
	double div_rate = 0;				// divergence rate of trajectories defined as log(d/d0).
	double le=0;						// Lyapunov exponent defined as the average of div_rate along the trajectories.
	int T = 100000;
	// time_series_file.open("test_data_lypo_r45.d");
	for (int i = 1; i <= T; ++i)
	{
		Y = X+E;						// Initial condition for the second trajectory

		//------Simultaneous evolution of the two trajectories by one time step ------
		X = RK4_integrator(X, dt, F);
		Y = RK4_integrator(Y, dt, F);
		//---------------------------------------------------------------------------

		E = Y-X;						// The tangent vector after one time step
		// cout<<E<<endl;
		d = Euclidean_measure(E);		// Length of the tangent vector
		// cout<<"d: "<<d<<endl;
		div_rate = log(d/d0);			// Local divergence rate of the trajectories 
		le +=div_rate;					// For calculating running average of divergence rate	
		E = (d0/d)*E;					// Renormalization of the tangent vector to length d0
		// cout<<E<<endl;
		// time_series_file<<div_rate<<endl;	// Writing local div rate in file
		if (i%1000==0)
		{
			cout<<i<<"\t"<<le/i<<endl;		//  Writing running average on std output 
		}
	}
	// time_series_file.close();
	le = le/(T*dt);
	cout<<"Lyapunov exponent: \t"<<le<<endl;		//  Writing running average on std output 

	// =========================================================================================
	return le;

}
