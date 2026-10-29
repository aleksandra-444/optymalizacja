#include"user_funs.h"

matrix ff0T(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla przypadku testowego
{
	matrix y;												// y zawiera warto�� funkcji celu
	y = pow(x(0) - ud1(0), 2) + pow(x(1) - ud1(1), 2);		// ud1 zawiera wsp�rz�dne szukanego optimum
	return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2)				// funkcja celu dla problemu rzeczywistego
{
	matrix y;												// y zawiera warto�� funkcji celu
	matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz�tkowe
		MT = matrix(2, new double[2] { m2d(x), 0.5 });		// MT zawiera moment si�y dzia�aj�cy na wahad�o oraz czas dzia�ania
	matrix* Y = solve_ode(df0, 0, 0.1, 20, Y0, ud1, MT);	// rozwi�zujemy r�wnanie r�niczkowe
	int n = get_len(Y[0]);									// d�ugo�� rozwi�zania
	double teta_max = Y[1](0, 0);							// szukamy maksymalnego wychylenia wahad�a
	for (int i = 1; i < n; ++i)
		if (teta_max < Y[1](i, 0))
			teta_max = Y[1](i, 0);
	y = abs(teta_max - m2d(ud1));							// warto�� funkcji celu (ud1 to za�o�one maksymalne wychylenie)
	Y[0].~matrix();											// usuwamy z pami�ci rozwi�zanie RR
	Y[1].~matrix();
	return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2)
{
	matrix dY(2, 1);										// definiujemy wektor pochodnych szukanych funkcji
	double m = 1, l = 0.5, b = 0.5, g = 9.81;				// definiujemy parametry modelu
	double I = m * pow(l, 2);
	dY(0) = Y(1);																// pochodna z po�o�enia to pr�dko��
	dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;	// pochodna z pr�dko�ci to przyspieszenie
	return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2)				
{
	matrix y;												
	y = -cos(0.1*x())*exp(-pow(0.1*x()-2*M_PI,2))+0.002*pow(0.1*x,2);		
	return y;
}

matrix ff1R(matrix x, matrix ud1, matrix ud2)				
{
	matrix y;												
	matrix Y0 = matrix(3, new double[3] {5, 1, 20});								
			
	matrix* Y = solve_ode(df1, 0, 1, 2000, Y0, ud1, x);	
	int n = get_len(Y[0]);									
	double max_temp = Y[1](0, 2);							
	for (int i = 1; i < n; ++i)
		if (max_temp < Y[1](i, 2))
			max_temp = Y[1](i, 2);
	y = abs(max_temp - 50);							
	delete[]Y;
	return y;
}

matrix df1(double t, matrix Y, matrix ud1, matrix ud2)
{
	double a = 0.98, b = 0.63, g = 9.81;				//a - wspolczynnik lepkosci cieczy, b - wspolczynnik zwezenia strumienia cieczy, g - przyspieszenie ziemskie
	double PA = 1.5, TA = 90, PB = 1;					//PA-pole podstawy, TA - temperatura, PB - pole podstawy zbiornik B
	double DB = 0.00365665, Fin = 0.01, Tin = 20;		//DB - pole przekroju, Fin - szybkosc wlewania si� wody do zbiornika B, Tin - temperatura wody

	matrix dY(3, 1);

	double FAout = Y(0) > 0 ?
		a * b * m2d(ud2) * 1e-4 * sqrt(2 * g * Y(0) / PA) : 0;

	double FBout = Y(1) > 0 ?
		a * b * DB * sqrt(2 * g * Y(1) / PB) : 0;
	
	dY(0) = -FAout;																
	dY(1) = FAout + Fin - FBout;
	dY(2)= Fin/Y(1)*(Tin-Y(2))+FAout/Y(1)*(TA-Y(2));
	return dY;
}



matrix ff2T(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	y = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
	return y;
}


matrix ff2R(matrix x, matrix ud1, matrix ud2)
{
	matrix y;
	matrix Y0(2, 1), Y_ref(2, new double[2] {3.14, 0});

	matrix* Y = solve_ode(df2, 0, 0.1, 100, Y0, Y_ref, x);
	int n = get_len(Y[0]);
	y = 0;
	
	for (int i = 1; i < n; ++i) {
		y = y + 10 * pow(Y_ref(0) - Y[1](i, 0), 2) + pow(Y_ref(1) - Y[1](i, 1), 2) + pow(x(0) * (Y_ref(0) - Y[1](i, 0)) + x(1) * (Y_ref(1) - Y[1](i, 1)), 2);
	}
	y = y * 0.1;
	delete[] Y;
	
	return y;
}

matrix df2(double t, matrix Y, matrix ud1, matrix ud2)
{

	//
	double mr = 1, mc = 5, l = 2, b = 0.25;
	double I = mr * l * l / 3 + mc * l * l;
	matrix dY(2, 1);
	dY(0) = Y(1);
	dY(1) = (ud2(0) * (ud1(0) - Y(0)) + ud2(1) * (ud1(1) - Y(1)) - b * Y(1)) / I;
	return dY;
}