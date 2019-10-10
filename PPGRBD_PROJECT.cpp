//NAIVNI ALGORITAM
#include <iostream>
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

vector<double> CramersRule(vector<vector<int>> points) {

	Matrix3i delta;
	delta << points[0][0], points[1][0], points[2][0],
			 points[0][1], points[1][1], points[2][1],
			 points[0][2], points[1][2], points[2][2];

	Matrix3i delta1;
	delta1 << points[3][0], points[1][0], points[2][0],
		      points[3][1], points[1][1], points[2][1],
		      points[3][2], points[1][2], points[2][2];

	Matrix3i delta2;
	delta2 << points[0][0], points[3][0], points[2][0],
		      points[0][1], points[3][1], points[2][1],
		      points[0][2], points[3][2], points[2][2];

	Matrix3i delta3;
	delta3 << points[0][0], points[1][0], points[3][0],
		      points[0][1], points[1][1], points[3][1],
		      points[0][2], points[1][2], points[3][2];

	double lambda1, lambda2, lambda3;
	vector<double> lambdas;

	lambda1 = delta1.determinant() * 1.0 / delta.determinant();
	lambdas.push_back(lambda1);

	lambda2 = delta2.determinant() * 1.0 / delta.determinant();
	lambdas.push_back(lambda2);

	lambda3 = delta3.determinant() * 1.0 / delta.determinant();
	lambdas.push_back(lambda3);

	return lambdas;
}

int main() {

	//ovaj deo ucitava oreginalne tacke(A, B, C, D) figure koje zelimo da projektujemo
	vector<vector<int>> original_points;

	//pomocni int za ucitavanje
	int x;

	//pomocni vektor za ucitavanje
	std::vector<int> tmp;

	//brojac za pojedinacne tacke, kad prodje 3 broja resetuje se
	int counter = 0;

	//brojac za petlju 
	int brojac = 0;

	//broj tacaka koje unosi korisnik
	int n;

	cout << "unesite broj tacaka" << endl;
	cin >> n;

	while (brojac < (n*3)) {

		cin >> x;
		brojac++;
		counter++;
		tmp.push_back(x);

		if (counter == 3)
		{
			original_points.push_back(tmp);
			tmp.clear();
			counter = 0;
		}
	}

	//D = lambda1*A + lambda2*B + lambda3*C - D je linearna kombinacija ostale 3 tacke
	vector<double> lambdas;
	
	//resavamo sistem uz pomoc Kramera i dobijamo lambde
	lambdas = CramersRule(original_points);

	cout << lambdas[0] << " " << lambdas[1] << " " << lambdas[2] << endl;

	//matrica prelaska iz kanonskog oblika(tacke su A0, B0, C0, D0) u nas cetvorougao ABCD: kolone su lambda1*A, lambda2*B i lambda3*C
	Matrix3d P;
	P << original_points[0][0] * lambdas[0], original_points[1][0] * lambdas[1], original_points[2][0] * lambdas[2],
		 original_points[0][1] * lambdas[0], original_points[1][1] * lambdas[1], original_points[2][1] * lambdas[2],
		 original_points[0][2] * lambdas[0], original_points[1][2] * lambdas[1], original_points[2][2] * lambdas[2];

	cout << P << endl;

	//Provera da li nas P transformise iz D0 u D
	MatrixXd D0;
	D0.resize(3, 1);
	D0 << 1, 1, 1;
	MatrixXd D;

	D = P * D0;

	cout << D << endl;

	//Sve isto ko za oreginalne tacke - ovde se samo radi o tackama u koje zelimo da projektujemo nase oreginalne
	// tacke Ap, Bp, Cp, Dp
	vector<vector<int>> destination_points;
	int x2;
	std::vector<int> tmp2;
	int counter2 = 0, brojac2 = 0, n2;

	cout << "unesite broj tacaka po drugi put" << endl;
	cin >> n2;


	while (brojac2 < (n2*3)) {

		brojac2++;
		std::cin >> x2;
		counter2++;
		tmp2.push_back(x2);

		if (counter2 == 3)
		{
			destination_points.push_back(tmp2);
			tmp2.clear();
			counter2 = 0;
		}
	}

	vector<double> lambdas2;

	lambdas2 = CramersRule(destination_points);

	//matrica prelaska iz kanonskog oblika u tacke Ap, Bp, Cp, Dp
	Matrix3d Pp;
	Pp << destination_points[0][0] * lambdas2[0], destination_points[1][0] * lambdas2[1], destination_points[2][0] * lambdas2[2],
		  destination_points[0][1] * lambdas2[0], destination_points[1][1] * lambdas2[1], destination_points[2][1] * lambdas2[2],
		  destination_points[0][2] * lambdas2[0], destination_points[1][2] * lambdas2[1], destination_points[2][2] * lambdas2[2];

	cout << Pp << endl;

	//provera da li iz D0 preko P dolazimo do Dp
	MatrixXd Dp;

	Dp = Pp * D0;

	cout << Dp << endl;

	//matrica Pp*(P)^(-1) predstavlja matricu prelaska iz ABCD u ApBpCpDp
	Matrix3d Q;
	Q = Pp * P.inverse();

	cout << Q << endl;

	//Provera da li idemo iz D u Dp
	MatrixXd Provera;

	Provera = Q * D;

	cout << Provera << endl;

}