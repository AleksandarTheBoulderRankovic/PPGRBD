//NAIVNI ALGORITAM, DLT, INTERFEJS
#include <iostream>
//good(relativly) intro to eigen : https://www.youtube.com/watch?v=6mMjv-tA5Jk
#include <Eigen/Dense>
#include <vector>
#include <atlimage.h>
#include "CImg.h"
#include <Eigen/SVD>
#include <Eigen/Eigen>

using namespace cimg_library;
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


	//load an image

	//vector<vector<int>> click_points;

	//ovde cuvamo tacke na koje kliknemo u interfejsu
	vector<vector<int>> org_p;
	vector<int> click_tmp;

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
	for (vector<int> p : org_p) {
		for (int e : p) {
			std::cout << e << " ";
		}
		std::cout << endl;
	}
	
	//broj tacaka koje unosi korisnik
	int n;

	n = org_p.size();

	//sad imam unos tacaka preko interfejsa ali ako zatreba moze i obican unos
	/*
	//ovaj deo ucitava oreginalne tacke(A, B, C, D) figure koje zelimo da projektujemo
	//vector<vector<int>> org_p;
	
	//pomocni int za ucitavanje
	int x;

	//pomocni vektor za ucitavanje
	std::vector<int> tmp;

	//brojac za pojedinacne tacke, kad prodje 3 broja resetuje se
	int counter = 0;

	//brojac za petlju 
	int brojac = 0;

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

	//D = lambda1*A + lambda2*B + lambda3*C - D je linearna kombinacija ostale 3 tacke
	vector<double> lambdas;
	
	//resavamo sistem uz pomoc Kramera i dobijamo lambde
	lambdas = CramersRule(org_p);

	std::cout << lambdas[0] << " " << lambdas[1] << " " << lambdas[2] << endl;

	//matrica prelaska iz kanonskog oblika(tacke su A0, B0, C0, D0) u nas cetvorougao ABCD: kolone su lambda1*A, lambda2*B i lambda3*C
	Matrix3d P;
	P << org_p[0][0] * lambdas[0], org_p[1][0] * lambdas[1], org_p[2][0] * lambdas[2],
		 org_p[0][1] * lambdas[0], org_p[1][1] * lambdas[1], org_p[2][1] * lambdas[2],
		 org_p[0][2] * lambdas[0], org_p[1][2] * lambdas[1], org_p[2][2] * lambdas[2];

	std::cout << P << endl;

	//Provera da li nas P transformise iz D0 u D
	MatrixXd D0;
	D0.resize(3, 1);
	D0 << 1, 1, 1;
	MatrixXd D;

	D = P * D0;

	std::cout << D << endl;

	//Sve isto ko za oreginalne tacke - ovde se samo radi o tackama u koje zelimo da projektujemo nase oreginalne
	// tacke Ap, Bp, Cp, Dp
	vector<vector<int>> dst_p;
	int x2;
	std::vector<int> tmp2;
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


	for (vector<int> p : dst_p) {
		for (int e : p) {
			std::cout << e << " ";
		}
		std::cout << endl;
	}
	


	vector<double> lambdas2;

	lambdas2 = CramersRule(dst_p);

	//matrica prelaska iz kanonskog oblika u tacke Ap, Bp, Cp, Dp
	Matrix3d Pp;
	Pp << dst_p[0][0] * lambdas2[0], dst_p[1][0] * lambdas2[1], dst_p[2][0] * lambdas2[2],
		  dst_p[0][1] * lambdas2[0], dst_p[1][1] * lambdas2[1], dst_p[2][1] * lambdas2[2],
		  dst_p[0][2] * lambdas2[0], dst_p[1][2] * lambdas2[1], dst_p[2][2] * lambdas2[2];

	std::cout << Pp << endl;

	//provera da li iz D0 preko P dolazimo do Dp
	MatrixXd Dp;

	Dp = Pp * D0;

	std::cout << Dp << endl;

	//matrica Pp*(P)^(-1) predstavlja matricu prelaska iz ABCD u ApBpCpDp
	Matrix3d Q;
	Q = Pp * P.inverse();

	std::cout << Q << endl;

	//Provera da li idemo iz D u Dp
	MatrixXd Provera;

	Provera = Q * D;

	std::cout << Provera << endl;


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
	//ove tri matrice nam pomazu da se izgradi ona velika (broj_tacaka*2) X 9 matrica
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
			A.resize((i+1)*2, 9);
			A = C;
			C.resize((i+2)*2, 9);
			i++;
		
		}

		std::cout << A << endl;

		//SVD dekompozicija A = UDV.transpose() - nas zanima samo zadnja kolona matricе V
		BDCSVD<MatrixXd> svd(A, ComputeFullU | ComputeFullV);

		std::cout << svd.computeV() << endl;
		std::cout << svd.matrixV() << endl;

		MatrixXd V;
		Matrix3d DLT;

		V.resize(9, 9);

		V << svd.matrixV();

		//uzimamo zadnju kolonu matrice V delimo je sa V(0,0) i mnozimo sa Q(0, 0)
		DLT << V(0, 8) * Q(0, 0) / V(0, 8), V(1, 8) * Q(0, 0) / V(0, 8), V(2, 8) * Q(0, 0) / V(0, 8),
			V(3, 8) * Q(0, 0) / V(0, 8), V(4, 8) * Q(0, 0) / V(0, 8), V(5, 8) * Q(0, 0) / V(0, 8),
			V(6, 8) * Q(0, 0) / V(0, 8), V(7, 8) * Q(0, 0) / V(0, 8), V(8, 8) * Q(0, 0) / V(0, 8);

		//DLT matrica treba da bude slicna kao Q matrica koju smo dobili iz naivnog algoritma
		std::cout << DLT << endl;

		//Provera da li idemo iz D u Dp
		MatrixXd Provera3;

		Provera3 = DLT * D;

		std::cout << Provera3 << endl;

		return 0;
}