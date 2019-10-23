//NAIVNI ALGORITAM, DLT, INTERFEJS
#include <iostream>
//good(relativly) intro to eigen : https://www.youtube.com/watch?v=6mMjv-tA5Jk
#include <Eigen/Dense>
#include <vector>
#include <atlimage.h>
#include "CImg.h"
#include <Eigen/SVD>
#include <Eigen/Eigen>
#include <cmath>

using namespace cimg_library;
using namespace Eigen;
using namespace std;

vector<double> CramersRule(vector<vector<double>> &points) {



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

Matrix3d& DLT_ALG(vector<vector<double>> &org_p, vector<vector<double>> &dst_p, int n) {

	MatrixXd A;
	MatrixXd B;
	MatrixXd C;

	int i = 1;

	A.resize(2, 9);
	B.resize(2, 9);
	C.resize(4, 9);

	A << 0, 0, 0,
		-org_p[0][0] * dst_p[0][2], -org_p[0][1] * dst_p[0][2], -org_p[0][2] * dst_p[0][2],
		org_p[0][0] * dst_p[0][1], org_p[0][1] * dst_p[0][1], org_p[0][2] * dst_p[0][1],
		org_p[0][0] * dst_p[0][2], org_p[0][1] * dst_p[0][2], org_p[0][2] * dst_p[0][2],
		0, 0, 0,
		-org_p[0][0] * dst_p[0][0], -org_p[0][1] * dst_p[0][0], -org_p[0][2] * dst_p[0][0];

	while (i < n) {
		B << 0, 0, 0,
			-org_p[i][0] * dst_p[i][2], -org_p[i][1] * dst_p[i][2], -org_p[i][2] * dst_p[i][2],
			org_p[i][0] * dst_p[i][1], org_p[i][1] * dst_p[i][1], org_p[i][2] * dst_p[i][1],
			org_p[i][0] * dst_p[i][2], org_p[i][1] * dst_p[i][2], org_p[i][2] * dst_p[i][2],
			0, 0, 0,
			-org_p[i][0] * dst_p[i][0], -org_p[i][1] * dst_p[i][0], -org_p[i][2] * dst_p[i][0];

		C << A, B;
		A.resize((i + 1) * 2, 9);
		A = C;
		C.resize((i + 2) * 2, 9);
		i++;

	}
	cout << "ogromna matrica dimnzija (2*n, 9) na koju primenjujemo SVD: " << endl << endl;
	std::cout << A << endl << endl;

	//SVD dekompozicija A = UDV.transpose() - nas zanima samo zadnja kolona matricе V
	BDCSVD<MatrixXd> svd(A, ComputeFullU | ComputeFullV);

	//std::cout << svd.computeV() << endl;

	std::cout << "Matrica V dobijena SVD dekompozicijom: " << endl << endl;
	std::cout << svd.matrixV() << endl << endl;

	MatrixXd V;
	Matrix3d DLT;

	V.resize(9, 9);

	V << svd.matrixV();

	//uzimamo zadnju kolonu matrice V delimo je sa V(0,0) i mnozimo sa Q(0, 0)
	DLT << V(0, 8), V(1, 8), V(2, 8),
		V(3, 8), V(4, 8), V(5, 8),
		V(6, 8), V(7, 8), V(8, 8);

	return DLT;
}

Matrix3d& naivni_ALG(vector<vector<double>> &org_p, vector<vector<double>> &dst_p) {

	//D = lambda1*A + lambda2*B + lambda3*C - D je linearna kombinacija ostale 3 tacke
	vector<double> lambdas;

	//resavamo sistem uz pomoc Kramera i dobijamo lambde
	lambdas = CramersRule(org_p);

	cout << "lambde za oreginalne tacke: " << endl;
	cout << lambdas[0] << " " << lambdas[1] << " " << lambdas[2] << endl << endl;

	//matrica prelaska iz kanonskog oblika(tacke su A0, B0, C0, D0) u nas cetvorougao ABCD: kolone su lambda1*A, lambda2*B i lambda3*C
	Matrix3d P;
	P << org_p[0][0] * lambdas[0], org_p[1][0] * lambdas[1], org_p[2][0] * lambdas[2],
		org_p[0][1] * lambdas[0], org_p[1][1] * lambdas[1], org_p[2][1] * lambdas[2],
		org_p[0][2] * lambdas[0], org_p[1][2] * lambdas[1], org_p[2][2] * lambdas[2];

	cout << "matrica prelaska iz kanonskog oblika u oreginalne tacke:" << endl << endl;
	cout << P << endl << endl;


	//Provera da li nas P transformise iz D0 u D
	MatrixXd D0;
	D0.resize(3, 1);
	D0 << 1, 1, 1;
	MatrixXd D;

	D = P * D0;

	cout << "Provera da li nas P transformise iz D0 u D: " << endl;
	cout << D << endl << endl;

	vector<double> lambdas2;

	lambdas2 = CramersRule(dst_p);

	cout << "lambde za tacke slike: " << endl;
	cout << lambdas2[0] << " " << lambdas2[1] << " " << lambdas2[2] << endl << endl;

	//matrica prelaska iz kanonskog oblika u tacke Ap, Bp, Cp, Dp
	Matrix3d Pp;
	Pp << dst_p[0][0] * lambdas2[0], dst_p[1][0] * lambdas2[1], dst_p[2][0] * lambdas2[2],
		dst_p[0][1] * lambdas2[0], dst_p[1][1] * lambdas2[1], dst_p[2][1] * lambdas2[2],
		dst_p[0][2] * lambdas2[0], dst_p[1][2] * lambdas2[1], dst_p[2][2] * lambdas2[2];

	cout << "matrica prelaska iz kanonskog oblika u tacke slike:" << endl << endl;
	std::cout << Pp << endl << endl;

	//provera da li iz D0 preko P dolazimo do Dp
	MatrixXd Dp;

	Dp = Pp * D0;

	cout << "Provera da li nas Pp transformise iz D0 u Dp: " << endl;
	std::cout << Dp << endl << endl;

	//matrica Pp*(P)^(-1) predstavlja matricu prelaska iz ABCD u ApBpCpDp
	Matrix3d Q;
	Q = Pp * P.inverse();

	return Q;
}

int main() {


	//load an image

	//vector<vector<int>> click_points;

	//ovde cuvamo tacke na koje kliknemo u interfejsu
	vector<vector<double>> org_p;
	vector<double> click_tmp;

	//kreira se slika koju ucitavamo iz bmp fajla
	CImg<unsigned char> image("building.bmp");
	//pomocna slika zato sto na prvoj imam neko obelezavanje izabranih tacaka koje ne zelim da imam na konstruisanoj slici
	CImg<unsigned char> image2("building.bmp");

	//svaku sliku koja nam dodje reskaliramo da izgleda interfejs lepse
	int width = 600;
	int height = 400;
	image.resize(width, height);
	image2.resize(width, height);


	//pravimo novu sliku koja je iste visine i sirine kao oreginalna ona je skroz crna
	CImg<unsigned char> visu(width, height, 1, 3, 0);

	//pravimo prozore i ucitavamo slike u njih
	CImgDisplay original_picture(image, "Click 4 or more points - press Enter to continue"), draw_disp(visu, "The processed image");

	//pomeramo prozore da interfejs izgleda lepse
	original_picture.move(0, 30);
	draw_disp.move(width, 30);

	//dok ne pritisnemo enter ili zatvorimo prozore radi se unos tacaka (samo klikni negde na unetoj slici)
	while (!original_picture.is_keyENTER() && !draw_disp.is_closed() && !original_picture.is_closed()) {
		original_picture.wait();
		if (original_picture.button() && original_picture.mouse_y() >= 0) {

			//unos koordinata kliknute tacke
			const int y = original_picture.mouse_y();
			const int x = original_picture.mouse_x();

			//unosimo koordinate kliknutog piksela u privremeni vector
			click_tmp.push_back(x);
			click_tmp.push_back(y);
			click_tmp.push_back(1);
			//click_points.push_back(click_tmp);

			//ovde taj privremeni vektor upisujemo u vector oreginalnih_tacaka
			org_p.push_back(click_tmp);
			click_tmp.clear();


			std::cout << "(" << x << "," << y << ") ="
				<< " R" << (int)image(x, y, 0, 0)
				<< " G" << (int)image(x, y, 0, 1)
				<< " B" << (int)image(x, y, 0, 2) << endl;

			visu(x, y, 0, 0) = (int)image(x, y, 0, 0); //R
			visu(x, y, 0, 1) = (int)image(x, y, 0, 1); //G
			visu(x, y, 0, 2) = (int)image(x, y, 0, 2); //B

			for (int i = 0; i < 8; i++)
			{
				for (int j = 0; j < 8; j++) {
					image(x + i, y + j, 0, 0) = 0;//(int)image(x, y, 0, 0); //R
					image(x + i, y + j, 0, 1) = 0;// (int)image(x, y, 0, 1); //G
					image(x + i, y + j, 0, 2) = 255;// (int)image(x, y, 0, 2); //B

					image(x - i, y - j, 0, 0) = (int)image(x, y, 0, 0); //R
					image(x - i, y - j, 0, 1) = (int)image(x, y, 0, 1); //G
					image(x - i, y - j, 0, 2) = (int)image(x, y, 0, 2); //B

					image(x + i, y - j, 0, 0) = (int)image(x, y, 0, 0); //R
					image(x + i, y - j, 0, 1) = (int)image(x, y, 0, 1); //G
					image(x + i, y - j, 0, 2) = (int)image(x, y, 0, 2); //B

					image(x - i, y + j, 0, 0) = (int)image(x, y, 0, 0); //R
					image(x - i, y + j, 0, 1) = (int)image(x, y, 0, 1); //G
					image(x - i, y + j, 0, 2) = (int)image(x, y, 0, 2); //B
				}
			}
			image.display(original_picture);
			visu.display(draw_disp);

			}
	}

	cout << "Oreginalne tacke:" << endl;
	for (vector<double> p : org_p) {
		for (double e : p) {
			std::cout << e << " ";
		}
		std::cout << endl;
	}
	cout << endl;
	
	//broj tacaka koje unosi korisnik
	int n;

	n = org_p.size();
	
	//sad imam unos tacaka preko interfejsa ali ako zatreba moze i obican unos
	/*
	//ovaj deo ucitava oreginalne tacke(A, B, C, D) figure koje zelimo da projektujemo
	vector<vector<double>> org_p;
	
	//pomocni int za ucitavanje
	double x;

	//pomocni vektor za ucitavanje
	std::vector<double> tmp;

	//brojac za pojedinacne tacke, kad prodje 3 broja resetuje se
	int counter = 0;

	//brojac za petlju 
	int brojac = 0;
	int n;
	cout << "unesite broj tacaka" << endl;
	cin >> n;
	cout << "unesite oreginalne tacke: " << endl;
	while (brojac < (n*2)) {

		cin >> x;
		brojac++;
		counter++;
		tmp.push_back(x);

		if (counter == 2)
		{
			tmp.push_back(1);
			org_p.push_back(tmp);
			tmp.clear();
			counter = 0;
		}
	}
	*/
	//Sve isto ko za oreginalne tacke - ovde se samo radi o tackama u koje zelimo da projektujemo nase oreginalne
	// tacke Ap, Bp, Cp, Dp
	vector<vector<double>> dst_p;
	double x2;
	std::vector<double> tmp2;
	int counter2 = 0, brojac2 = 0;

	std::cout << "Unesite tacke slike: " << endl;


	while (brojac2 < (n*2)) {

		brojac2++;
		std::cin >> x2;
		counter2++;
		tmp2.push_back(x2);

		if (counter2 == 2)
		{
			tmp2.push_back(1);
			dst_p.push_back(tmp2);
			tmp2.clear();
			counter2 = 0;
		}
	}
	
	Matrix3d Q = naivni_ALG(org_p, dst_p);

	cout << "matrica transformacije dobijena naivnim algoritmom: " << endl << endl;
	cout << Q << endl << endl;
	
	//Provera da li idemo iz D u Dp
	MatrixXd Provera;
	MatrixXd D;
	D.resize(3, 1);
	D << org_p[3][0], org_p[3][1], org_p[3][2];

	Provera = Q * D;

	std::cout << "Provera da li idemo iz D u Dp naivnim algoritmom" << endl << endl << Provera << endl << endl;
	
	//iscrtavanje izmenjene oreginalne slike
	MatrixXd TmpOrgPoint;
	TmpOrgPoint.resize(3, 1);

	MatrixXd TmpDstPoint;
	TmpDstPoint.resize(3, 1);

	int cp, rp, cpp, rpp;
	
	//dva nacina za konstrukciju izmenjene slike - kada iteriramo po pikselima oreginalne slike i koristimo Q matricu
	//a mozemo i da iteriramo kroz sliku i koristimo matricu Q.inverse()
	/*

	//prvi nacin za konstrukciju slike - daje malo losije rezultate
	for (int r = 0; r < height; r++)
		for (int c = 0; c < width; c++)
		{
			TmpOrgPoint << c, r, 1;

			TmpDstPoint = Q * TmpOrgPoint;

			if (TmpDstPoint(2, 0) != 0) {
				cp = floor(TmpDstPoint(0, 0) / TmpDstPoint(2, 0));
				rp = floor(TmpDstPoint(1, 0) / TmpDstPoint(2, 0));

				cpp = ceil(TmpDstPoint(0, 0) / TmpDstPoint(2, 0));
				rpp = ceil(TmpDstPoint(1, 0) / TmpDstPoint(2, 0));

			}
			else {
				continue;
			}
			if (rp > height || rp < 0 || cp > width || cp < 0)
				continue;
			if (rpp > height || rpp < 0 || cpp > width || cpp < 0)
				continue;
			
			visu(cp, rp, 0, 0) = (int)image2(c, r, 0, 0); //R
			visu(cp, rp, 0, 1) = (int)image2(c, r, 0, 1); //G
			visu(cp, rp, 0, 2) = (int)image2(c, r, 0, 2); //B

			visu(cpp, rpp, 0, 0) = (int)image2(c, r, 0, 0); //R
			visu(cpp, rpp, 0, 1) = (int)image2(c, r, 0, 1); //G
			visu(cpp, rpp, 0, 2) = (int)image2(c, r, 0, 2); //B
		}
	*/

	//invertovana Q matrica
	Matrix3d Qinv;
	Qinv = Q.inverse();
	
	//drugi nacin za konstrukciju slike
	//iteriramo kroz sve piksele slike koju konstruisemo
	for (int r = 0; r < height; r++)
		for (int c = 0; c < width; c++)
		{
			TmpDstPoint << c, r, 1;

			//trzimo koordinate pikela oreginalne slike koji odgovaraju trenutnom pikselu slike koju konstruisemo
			TmpOrgPoint = Qinv * TmpDstPoint;

			//pretvaramo homogene koordinate piksela u afine
			if (TmpOrgPoint(2, 0) != 0) {
				cp = floor(TmpOrgPoint(0, 0) / TmpOrgPoint(2, 0));
				rp = floor(TmpOrgPoint(1, 0) / TmpOrgPoint(2, 0));

				cpp = ceil(TmpOrgPoint(0, 0) / TmpOrgPoint(2, 0));
				rpp = ceil(TmpOrgPoint(1, 0) / TmpOrgPoint(2, 0));
			}
			else {
				continue;
			}
			if (rp > height || rp < 0 || cp > width || cp < 0)
				continue;
			if (rpp > height || rpp < 0 || cpp > width || cpp < 0)
				continue;

			//iscrtavamo piksel na odgovarajucoj lokaciji sa istom bojom kao u oreginalnoj slici
			visu(c, r, 0, 0) = (int)image2(cp, rp, 0, 0); //R
			visu(c, r, 0, 1) = (int)image2(cp, rp, 0, 1); //G
			visu(c, r, 0, 2) = (int)image2(cp, rp, 0, 2); //B

			visu(c, r, 0, 0) = (int)image2(cpp, rpp, 0, 0); //R
			visu(c, r, 0, 1) = (int)image2(cpp, rpp, 0, 1); //G`	
			visu(c, r, 0, 2) = (int)image2(cpp, rpp, 0, 2); //B

		}

	visu.display(draw_disp);
	while (!draw_disp.is_closed())
		draw_disp.wait();
	
	//DLT

		//DLT matrica treba da bude slicna kao Q matrica koju smo dobili iz naivnog algoritma
		Matrix3d V = DLT_ALG(org_p, dst_p, n);
		cout << "Matrica dobijena DLT algoritmom: " << endl << endl;
		cout << V << endl << endl;

		Matrix3d DLT_skalirano;
		//uzimamo zadnju kolonu matrice V delimo je sa V(0,0) i mnozimo sa Q(0, 0)
		DLT_skalirano << V(0, 0)* Q(0, 0) / V(0, 0), V(0, 1)* Q(0, 0) / V(0, 0), V(0, 2)* Q(0, 0) / V(0, 0),
			V(1, 0)* Q(0, 0) / V(0, 0), V(1, 1)* Q(0, 0) / V(0, 0), V(1, 2)* Q(0, 0) / V(0, 0),
			V(2, 0)* Q(0, 0) / V(0, 0), V(2, 1)* Q(0, 0) / V(0, 0), V(2, 2)* Q(0, 0) / V(0, 0);

		cout << "Matrica transformacije dobijena DLT-om ali skalirana:" << endl << endl;
		cout << DLT_skalirano << endl << endl;


		//Provera da li idemo iz D u Dp
		MatrixXd Provera3;

		Provera3 = DLT_skalirano * D;

		cout << "provera da li idemo iz D u Dp sa skaliranom DLT matricom:" << endl << endl;
		cout << Provera3 << endl;

	//Normalizovani DLT

		//Prvo trazimo teziste sistema tecaka (oreginalnih i slike)
		int sum_x = 0, sum_y = 0, sum_x_s = 0, sum_y_s = 0;
		for (int i = 0; i < n; i++) {
			sum_x += org_p[i][0];

			sum_y += org_p[i][1];

			sum_x_s += dst_p[i][0];

			sum_y_s += dst_p[i][1];
		}

		vector<double> centar;
		centar.push_back(sum_x*1.0 / n);
		centar.push_back(sum_y*1.0 / n);

		vector<double> centar_s;
		centar_s.push_back(sum_x_s * 1.0 / n);
		centar_s.push_back(sum_y_s * 1.0 / n);


		//koordinate centra za oreginal i sliku:
		cout << "centri oreginala i slike" << endl;
		cout << centar[0] << " " << centar[1] << endl;

		cout << centar_s[0] << " " << centar_s[1] << endl << endl;

		//racunamo prosecnu razdaljinu drugih tacaka od centra
		double sum_len=0, sum_len_s = 0;
		
		for (int i = 0; i < n; i++) {
			sum_len += sqrt((org_p[i][0]-centar[0])* (org_p[i][0] - centar[0]) + (org_p[i][1] - centar[1])*(org_p[i][1] - centar[1]));

			sum_len_s += sqrt((dst_p[i][0] - centar_s[0]) * (dst_p[i][0] - centar_s[0]) + (dst_p[i][1] - centar_s[1]) * (dst_p[i][1] - centar_s[1]));
			
		}
		cout << "prosecno rastojanje od centra kod oreginalnih tacaka:" << endl;
		sum_len = sum_len / n;
		cout << sum_len << endl << endl;
		

		cout << "prosecno rastojanje od centra kod tacaka slike:" << endl;
		sum_len_s = sum_len_s / n;
		cout << sum_len_s << endl << endl;

		//sad pravimo matrice T i Tp
		Matrix3d T;

		T << (sqrt(2) / sum_len), 0, -centar[0],
			 0, (sqrt(2) / sum_len), -centar[1],
			 0, 0, 1;

		Matrix3d Tp;

		Tp << (sqrt(2) / sum_len_s), 0, -centar_s[0],
			0, (sqrt(2) / sum_len_s), -centar_s[1],
			0, 0, 1;

		//transformacija koordinata oreginala i slika u normalizovane koordinate
		cout << "Matrice T i Tp: " << endl << T << endl << endl << Tp << endl << endl;
		MatrixXd tmpPoint;
		MatrixXd tmpNormPoint;
		tmpPoint.resize(3, 1);

		vector<vector <double>> norm_org_p;
		vector<double> tmpVec;
		for (int i = 0; i < n; i++) {
			tmpPoint << org_p[i][0], org_p[i][1], org_p[i][2];
			tmpNormPoint = T * tmpPoint;
			tmpVec.push_back(tmpNormPoint(0, 0));
			tmpVec.push_back(tmpNormPoint(1, 0));
			tmpVec.push_back(tmpNormPoint(2, 0));
			norm_org_p.push_back(tmpVec);
			tmpVec.clear();	
		}

		cout << "nove normalizovane tacke oreginala:" << endl << endl;
		for (vector<double> aca : norm_org_p) {
			for (double a : aca) {
				cout << a << " ";
			}
			cout << endl;
		}
		cout << endl;

		vector<vector <double>> norm_org_p_s;

		for (int i = 0; i < n; i++) {
			tmpPoint << dst_p[i][0], dst_p[i][1], dst_p[i][2];
			tmpNormPoint = Tp * tmpPoint;
			tmpVec.push_back(tmpNormPoint(0, 0));
			tmpVec.push_back(tmpNormPoint(1, 0));
			tmpVec.push_back(tmpNormPoint(2, 0));
			norm_org_p_s.push_back(tmpVec);
			tmpVec.clear();
		}
		cout << "nove normalizovane tacke slike:" << endl << endl;
		for (vector<double> aca : norm_org_p_s) {
			for (double a : aca) {
				cout << a << " ";
			}
			cout << endl;
		}
		cout << endl;
		
		Matrix3d normDLT = DLT_ALG(norm_org_p, norm_org_p_s, n);

		cout << "Matrica dobijena primenom DLT-a nad normiranim tackama:" << endl << normDLT << endl << endl;

		Matrix3d lastMatrix;

		lastMatrix = Tp.inverse() * normDLT * T;

		cout << "matrica dobijena normalizovanim DLT algoritmom:" << endl << lastMatrix << endl << endl;

		Matrix3d skalirano_normirano_DLT;

		skalirano_normirano_DLT << lastMatrix(0, 0)* Q(0, 0) / lastMatrix(0, 0)
			, lastMatrix(0, 1)* Q(0, 0) / lastMatrix(0, 0)
			, lastMatrix(0, 2)* Q(0, 0) / lastMatrix(0, 0)
			, lastMatrix(1, 0)* Q(0, 0) / lastMatrix(0, 0)
			, lastMatrix(1, 1)* Q(0, 0) / lastMatrix(0, 0)
			, lastMatrix(1, 2)* Q(0, 0) / lastMatrix(0, 0)
			, lastMatrix(2, 0)* Q(0, 0) / lastMatrix(0, 0)
			, lastMatrix(2, 1)* Q(0, 0) / lastMatrix(0, 0)
			, lastMatrix(2, 2)* Q(0, 0) / lastMatrix(0, 0);

		cout << "poslednja stvar ikada:" << endl << endl << skalirano_normirano_DLT << endl << endl;


		return 0;
}