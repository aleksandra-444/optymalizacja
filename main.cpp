/*********************************************
Kod stanowi uzupe³nienie materia³ów do æwiczeñ
w ramach przedmiotu metody optymalizacji.
Kod udostêpniony na licencji CC BY-SA 3.0
Autor: dr in¿. £ukasz Sztangret
Katedra Informatyki Stosowanej i Modelowania
Akademia Górniczo-Hutnicza
Data ostatniej modyfikacji: 15.09.2026
*********************************************/

#include"opt_alg.h"

#include <random>
#include <fstream>
#include <string>

void lab0();
void lab1();
void lab1_test();
void lab2();
void lab3();
void lab4();
void lab5();
void lab6();

int main()
{
	try
	{
		// lab1();  // Pelne laboratorium
		//lab1_test(); // Test kontrolny zbiornikow
        lab2();
    }
	catch (string EX_INFO)
	{
		cerr << "ERROR:\n";
		cerr << EX_INFO << endl << endl;
	}
	return 0;
}

void lab0()
{
	double epsilon;												// dok³adnoœæ
	int Nmax;													// maksymalna liczba wywo³añ funkcji celu
	matrix lb, ub;												// dolne oraz górne ograniczenie
	solution opt;												// rozwi¹zanie optymalne znalezione przez algorytm
	if (0) // Funkcja testowa
	{
		epsilon = 1e-2;
		Nmax = 10000;
		lb = matrix(2, 1, -5);
		ub = matrix(2, 1, 5);
		matrix a(2, 1);											// dok³adne rozwi¹zanie optymalne
		a(0) = -1;
		a(1) = 2;
		opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);			// wywo³anie procedury optymalizacji
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();								// wyzerowanie liczników
		int* n = get_size(opt.ud);
		cout << n[0] << '\t' << n[1] << endl << endl;
	}
	if (1) // Wahad³o
	{
		Nmax = 1000;
		epsilon = 1e-2;
		lb = 0, ub = 5;
		double teta_opt = 1;									// maksymalne wychylenie wahad³a
		opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);		// wywo³anie procedury optymalizacji
		cout << opt << endl << endl;							// wypisanie wyniku
		solution::clear_calls();								// wyzerowanie liczników

		// Zapis zmian wartoœci funckji celu do pliku
		ofstream Sout("optymalicja_lab0.csv");					// definiujemy strumieñ do pliku .csv
		Sout << opt.ud;											// zapisujemy wyniki w pliku
		Sout.close();											// zamykamy strumieñ

		// Zapis symulacji do pliku csv
		matrix Y0 = matrix(2, 1),								// Y0 zawiera warunki pocz¹tkowe
			MT = matrix(2, new double[2] { m2d(opt.x), 0.5 });	// MT zawiera moment si³y dzia³aj¹cy na wahad³o oraz czas dzia³ania
		matrix* Y = solve_ode(df0, 0, 0.1, 20, Y0, NAN, MT);	// rozwi¹zujemy równanie ró¿niczkowe
		Sout.open("symulacja_lab0.csv");						// otwieramy plik .csv
		Sout << hcat(Y[0], Y[1]);								// zapisyjemy wyniki w pliku
		Sout.close();											// zamykamy strumieñ
		Y[0].~matrix();											// usuwamy z pamiêci rozwi¹zanie RR
		Y[1].~matrix();
	}
}


void lab1_test()
{
    const double DA_cm2 = 50.0;
    const double temperatura_oczekiwana = 61.2859;

    matrix Y0(3, new double[3]{5.0, 1.0, 20.0});
    matrix* Y = solve_ode(df1, 0.0, 1.0, 2000.0, Y0, NAN, matrix(DA_cm2));

    int n = get_len(Y[0]);
    double Tmax = Y[1](0, 2);
    double tmax = Y[0](0);

    for (int i = 1; i < n; ++i)
    {
        if (Y[1](i, 2) > Tmax)
        {
            Tmax = Y[1](i, 2);
            tmax = Y[0](i);
        }
    }

    cout << "Test kontrolny: DA = " << DA_cm2 << " cm2\n";
    cout << "Maksymalna temperatura B: " << Tmax << " C\n";
    cout << "Czas wystapienia maksimum: " << tmax << " s\n";
    cout << "Wartosc oczekiwana: " << temperatura_oczekiwana << " C\n";
    cout << "Roznica: " << abs(Tmax - temperatura_oczekiwana) << " C\n";

    delete[] Y;
}

void lab1()
{
    double epsilon = 1e-2, gamma = 1e-4;
    int Nmax = 10000;
    double alpha[3] = { 1.8, 2.7, 4.2 };

    mt19937 gen(random_device{}());
    uniform_real_distribution<double> losuj(-100, 100);

    ofstream Sout("wyniki_lab1.csv");

    Sout << "alpha;proba;x0;a;b;f_exp;"
        << "x_fib;y_fib;f_fib;rodzaj_fib;"
        << "x_lag;y_lag;f_lag;rodzaj_lag\n";

    for (int j = 0; j < 3; ++j)
    {
        for (int i = 0; i < 100; ++i)
        {
            double x0 = losuj(gen);
            double a = 0, b = 0;
            int f_exp = 0;

            solution opt_fib, opt_lag;
            int f_fib = 0, f_lag = 0;

            string rodzaj_fib = "Brak zbieznosci";
            string rodzaj_lag = "Brak zbieznosci";

            // Metoda ekspansji
            try
            {
                solution::clear_calls();

                double* przedzial = expansion(
                    ff1T, x0, 1, alpha[j], Nmax
                );

                a = przedzial[0];
                b = przedzial[1];

                delete[] przedzial;

                f_exp = solution::f_calls;
            }
            catch (string ex_info)
            {
                Sout << alpha[j] << ';' << i + 1 << ';'
                    << x0 << ";;;0;;;;Brak zbieznosci"
                    << ";;;;Brak zbieznosci\n";
                continue;
            }

            // Metoda Fibonacciego
            try
            {
                solution::clear_calls();

                opt_fib = fib(ff1T, a, b, epsilon);
                f_fib = solution::f_calls;

                if (opt_fib.flag == 1)
                {
                    if (m2d(opt_fib.y) < -0.1)
                        rodzaj_fib = "Globalne";
                    else
                        rodzaj_fib = "Lokalne";
                }
            }
            catch (string ex_info)
            {
                f_fib = solution::f_calls;
            }

            // Metoda Lagrange'a
            try
            {
                solution::clear_calls();

                opt_lag = lag(
                    ff1T, a, b, epsilon, gamma, Nmax
                );

                f_lag = solution::f_calls;

                if (opt_lag.flag == 1)
                {
                    if (m2d(opt_lag.y) < -0.1)
                        rodzaj_lag = "Globalne";
                    else
                        rodzaj_lag = "Lokalne";
                }
            }
            catch (string ex_info)
            {
                f_lag = solution::f_calls;
            }

            // Zapis wynikow do Tabeli 1
            Sout << alpha[j] << ';'
                << i + 1 << ';'
                << x0 << ';'
                << a << ';'
                << b << ';'
                << f_exp << ';';

            if (rodzaj_fib != "Brak zbieznosci")
                Sout << m2d(opt_fib.x) << ';'
                << m2d(opt_fib.y);
            else
                Sout << ';';

            Sout << ';' << f_fib << ';'
                << rodzaj_fib << ';';

            if (rodzaj_lag != "Brak zbieznosci")
                Sout << m2d(opt_lag.x) << ';'
                << m2d(opt_lag.y);
            else
                Sout << ';';

            Sout << ';' << f_lag << ';'
                << rodzaj_lag << '\n';
        }
    }

    Sout.close();

    // Optymalizacja bez ekspansji
    ofstream Bez("bez_ekspansji_lab1.csv");
    Bez << "metoda;x;y;f_calls;flag\n";

    solution::clear_calls();
    solution opt = fib(ff1T, -100, 100, epsilon);

    Bez << "Fibonacci;"
        << m2d(opt.x) << ';'
        << m2d(opt.y) << ';'
        << solution::f_calls << ';'
        << opt.flag << '\n';

    ofstream Historia("przedzial_fib_lab1.csv");
    Historia << "iteracja;dlugosc\n";

    for (int i = 0; i < get_len(opt.ud); ++i)
        Historia << i + 1 << ';' << opt.ud(i) << '\n';

    Historia.close();

    solution::clear_calls();
    opt = lag(ff1T, -100, 100, epsilon, gamma, Nmax);

    Bez << "Lagrange;"
        << m2d(opt.x) << ';'
        << m2d(opt.y) << ';'
        << solution::f_calls << ';'
        << opt.flag << '\n';

    Historia.open("przedzial_lag_lab1.csv");
    Historia << "iteracja;dlugosc\n";

    for (int i = 0; i < get_len(opt.ud); ++i)
        Historia << i + 1 << ';' << opt.ud(i) << '\n';

    Historia.close();
    Bez.close();

    // Problem rzeczywisty - zbiorniki
    ofstream Real("wyniki_zbiorniki_lab1.csv");
    Real << "metoda;DA_cm2;y;f_calls;flag\n";

    solution::clear_calls();
    solution opt_fib = fib(ff1R, 1, 100, epsilon);
    int f_fib = solution::f_calls;

    Real << "Fibonacci;"
        << m2d(opt_fib.x) << ';'
        << m2d(opt_fib.y) << ';'
        << f_fib << ';'
        << opt_fib.flag << '\n';

    solution::clear_calls();
    solution opt_lag = lag(
        ff1R, 1, 100, epsilon, gamma, Nmax
    );
    int f_lag = solution::f_calls;

    Real << "Lagrange;"
        << m2d(opt_lag.x) << ';'
        << m2d(opt_lag.y) << ';'
        << f_lag << ';'
        << opt_lag.flag << '\n';

    Real.close();

    // Symulacja dla obu metod
    matrix Y0(3, new double[3] { 5, 1, 20 });

    matrix* Yfib = solve_ode(
        df1, 0, 1, 2000, Y0, NAN, opt_fib.x
    );

    matrix* Ylag = solve_ode(
        df1, 0, 1, 2000, Y0, NAN, opt_lag.x
    );

    ofstream Sym("symulacja_lab1.csv");

    Sym << "t;VA_fib;VA_lag;VB_fib;VB_lag;TB_fib;TB_lag\n";

    int n = get_len(Yfib[0]);

    for (int i = 0; i < n; ++i)
    {
        Sym << Yfib[0](i) << ';'
            << Yfib[1](i, 0) << ';'
            << Ylag[1](i, 0) << ';'
            << Yfib[1](i, 1) << ';'
            << Ylag[1](i, 1) << ';'
            << Yfib[1](i, 2) << ';'
            << Ylag[1](i, 2) << '\n';
    }

    Sym.close();

    delete[] Yfib;
    delete[] Ylag;

    cout << "Laboratorium 1 zakonczone." << endl;
    cout << "Wyniki zapisano do plikow CSV." << endl;
}


void lab2()
{

    //deklaracja parametrow optymalizacji
    //s - poczatkowa dlugosc kroku
    //alphaHJ - wspolczynnik zmniejszania kroku w metodzie hooke'a jeevesa
    //alphaR - wspolczynnik zwiekszania kroku w metodzie rosenbrocka
    //beta - wspolczynnik zmniejszania i odwracania krokuw metodzie rosenbrocka
    //epsilon - dokladnosc zakonczenia optymalizacji
    double s, alphaHJ, alphaR, beta, epsilon;
    int Nmax;//maksymalna liczba wywolan funkcji celu 

    //obiekt przechowywujacy wynik optymalizacji
    //opt.x - znaleziony punkt, opt.y - wartosc funkcji celu
    solution opt;
    //x0 - wektor punktu poczatkowego
    //s0 - wektor poczatkowych krokow dla metdoy rosenbrocka
    matrix x0, s0;

    if (0)//funkcja testowa
    {

        s = 0.5;//poczatkowa dlugosc kroku
        alphaHJ = 0.5;//hooke-jeeves: po nieudanym poszukiwaniu krok zmniejsza sie
        alphaR = 2;//rosenbrock: po udanym ruchu krok zwieksza sie dwukrotnie
        beta = 0.5;//rosenbrock: po nieudanym ruchu krok zmniejsza sie o polowe i zmienia znak czyli kierunek
        epsilon = 1e-3;//dokladnosc optymalizacji 0.001
        Nmax = 1000;//limit wywolan funkcji celu
        s0 = matrix(2, 1, s);//utworzenie wektora 2x1 zawierajacego poczatkowe kroki s0=[0.5;0.5]

        //losowanie punktu poczatkowego x0=[x1,x2]
        //rand_mat(2,1) generuje dwie liczby z przedzialu [0,1]
        //po pomnozeniu przez 2 i odjeciu 1 otzrymujemy x1,x2 z [-1,1]
        x0 = 2 * rand_mat(2, 1) - 1;
        cout << x0 << endl << endl;

        //metoda  hooke'a jeevsa
        //ff2T - testowa funkcje celu dwoch zmiennych
        //x0 - punkt startowy, s - dlugosc kroku
        //alphaHJ - wspolczynnik redukcji kroku
        // epsilon - dokladnosc , Nmax - limit wywolan
        //wynik zostaje zapisany w obiekcie opt
        opt = HJ(ff2T, x0, s, alphaHJ, epsilon, Nmax);
        cout << opt << endl << endl;
        //wyzerowanie licznikow wywolan funkcji celu przed rozpoczeciem metody rosenbrocka
        solution::clear_calls();
        //metoda rosenbrocka dla tej samej funkcji celu i punktu x0
        //s0 - wektor krokow, alphaR - ekspansja, beta - kontrakcja 
        //wynik rosenbrocka zastepuje poprzednhi wynik w opt
        opt = Rosen(ff2T, x0, s0, alphaR,beta, epsilon, Nmax);
        cout << opt << endl << endl;
        //wyserowanie
        solution::clear_calls();


    }
    if (0)//ramie robota
    {
        s = 0.5;//poczatkowa dlugosc kroku dla optymalizacji regulatora
        alphaHJ = 0.5;//zmniejszenie kroku o polowe
        alphaR = 2;//podwajanie kroku po udanym ruchu
        beta = 0.5;//zmniejszanie i odwracanie kroku po nieudanym ruchu
        epsilon = 1e-3;//dokladnosc
        Nmax = 1000;//limit wywolan
        x0 = 20 * rand_mat(2, 1);//losowanie poczatkowych wspolczynnikowregulatora x0(0)=k1, x0(1)=k2 oba naleza do przedzialo [0,20]

        cout << x0 << endl << endl;
        //optymalizacja refukatora metoda hookea jeevesa
        //ff2r funkcja celu dla modelu ramienia robota
        //szukamy wspolczynnikow k1, k2 minimalizujacych
        //funkcje jakosci sterowania Q(k1,k2)
        opt = HJ(ff2T, x0, s, alphaHJ, epsilon, Nmax);
        cout << opt << endl << endl;
        solution::clear_calls();
        //utworzenie wektora poczatkowych krkow dla rosenbrocka
        //s0=[2,2]
        s0 = matrix(2, 1, s);
        //optymalizacja tego samego regulatora metoda Rosenbrocka
        //obie metody rozpoczynaja obliczenia od tego samego punktu x0
        //rosenbrock moze niezaleznie zzmieniac kroki w obu kierunkach
        opt = Rosen(ff2T, x0, s0, alphaR,beta, epsilon, Nmax);
        cout << opt << endl << endl;
        solution::clear_calls();

        
    }
    if (1)//test modelu i funkcji celu
    {

        //ustawienie wspolczynnikow regulatora
        //k1=5, k2=5
        //matrix(2,1,5) tworzy wektor [5;5]
        //nie wykonujemy tutaj optymalizacji, tylko test funkcji celu
        opt.x = matrix(2, 1, 5);
        cout << endl << endl;
        //obliczenie i wyswietlenie wartosci funkcji celu Q(5,5)
        //fit_fun(ff2R) wywoluje funkcje celu dla punktu opt.x
        //ff2R przeprowadza symulacje ruchu ramienia robota
        //i oblicza funkcje jakosci sterowania
        //oczekiwany wynik wedlug konspektu okolo 775.229
        cout << opt.fit_fun(ff2R) << endl << endl;
        solution::clear_calls();
    }
}

void lab3()
{

}

void lab4()
{

}

void lab5()
{

}

void lab6()
{

}
