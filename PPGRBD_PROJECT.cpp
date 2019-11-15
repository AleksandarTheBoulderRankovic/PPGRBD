//NAIVNI ALGORITAM, DLT, MODIFIKOVANI DLT, INTERFEJS
#include <iostream>
//good(relativly) intro to eigen : https://www.youtube.com/watch?v=6mMjv-tA5Jk
#include <Eigen/Dense>
#include <vector>
#include <atlimage.h>
#include "CImg.h"
#include <Eigen/SVD>
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

Matrix3d& naivni_ALG(vector<vector<double>>& org_p, vector<vector<double>>& dst_p) {

	//D = lambda1*A + lambda2*B + lambda3*C - D je linearna kombinacija ostale 3 tacke
	vector<double> lambdas;

	//resavamo sistem uz pomoc Kramera i dobijamo lambde
	lambdas = CramersRule(org_p);

	cout << "lambde za originalne tacke: " << endl;
	cout << lambdas[0] << " " << lambdas[1] << " " << lambdas[2] << endl << endl;

	//matrica prelaska iz kanonskog oblika(tacke su A0, B0, C0, D0) u nas cetvorougao ABCD: kolone su lambda1*A, lambda2*B i lambda3*C
	Matrix3d P;
	P << org_p[0][0] * lambdas[0], org_p[1][0] * lambdas[1], org_p[2][0] * lambdas[2],
		org_p[0][1] * lambdas[0], org_p[1][1] * lambdas[1], org_p[2][1] * lambdas[2],
		org_p[0][2] * lambdas[0], org_p[1][2] * lambdas[1], org_p[2][2] * lambdas[2];

	cout << "matrica prelaska iz kanonskog oblika u originalne tacke:" << endl << endl;
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


Matrix3d& DLT_ALG(vector<vector<double>> &org_p, vector<vector<double>> &dst_p, int n) {
	
	//pomocne matrice od kojih pravimo Veliku matricu za SVD
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

	std::cout << "Matrica V dobijena SVD dekompozicijom: " << endl << endl;
	std::cout << svd.matrixV() << endl << endl;

	MatrixXd V;
	Matrix3d DLT;

	V.resize(9, 9);

	V << svd.matrixV();

	//uzimamo zadnju kolonu matrice V ona je nasa matrica transformacije
	DLT << V(0, 8), V(1, 8), V(2, 8),
		V(3, 8), V(4, 8), V(5, 8),
		V(6, 8), V(7, 8), V(8, 8);

	return DLT;
}

Matrix3d& normalizovani_DLT(vector<vector<double>> &org_p, vector<vector<double>> &dst_p, int n)
{
	//Prvo trazimo teziste sistema tecaka (originalnih i slike)
	int sum_x = 0, sum_y = 0, sum_x_s = 0, sum_y_s = 0;
	for (int i = 0; i < n; i++) {
		sum_x += org_p[i][0];

		sum_y += org_p[i][1];

		sum_x_s += dst_p[i][0];

		sum_y_s += dst_p[i][1];
	}

	vector<double> centar;
	centar.push_back(sum_x * 1.0 / n);
	centar.push_back(sum_y * 1.0 / n);

	vector<double> centar_s;
	centar_s.push_back(sum_x_s * 1.0 / n);
	centar_s.push_back(sum_y_s * 1.0 / n);


	//koordinate centra za original i sliku:
	cout << "centri originala i slike" << endl;
	cout << centar[0] << " " << centar[1] << endl;

	cout << centar_s[0] << " " << centar_s[1] << endl << endl;

	//racunamo prosecnu razdaljinu drugih tacaka od centra
	double sum_len = 0, sum_len_s = 0;

	for (int i = 0; i < n; i++) {
		sum_len += sqrt((org_p[i][0] - centar[0]) * (org_p[i][0] - centar[0]) + (org_p[i][1] - centar[1]) * (org_p[i][1] - centar[1]));

		sum_len_s += sqrt((dst_p[i][0] - centar_s[0]) * (dst_p[i][0] - centar_s[0]) + (dst_p[i][1] - centar_s[1]) * (dst_p[i][1] - centar_s[1]));

	}
	cout << "prosecno rastojanje od centra kod originalnih tacaka:" << endl;
	sum_len = sum_len / n;
	cout << sum_len << endl << endl;


	cout << "prosecno rastojanje od centra kod tacaka slike:" << endl;
	sum_len_s = sum_len_s / n;
	cout << sum_len_s << endl << endl;

	//sad pravimo matrice T i Tp koje normalizuju tacke originala i slike
	Matrix3d T;

	T << (sqrt(2) / sum_len), 0, -centar[0],
		0, (sqrt(2) / sum_len), -centar[1],
		0, 0, 1;

	Matrix3d Tp;

	Tp << (sqrt(2) / sum_len_s), 0, -centar_s[0],
		0, (sqrt(2) / sum_len_s), -centar_s[1],
		0, 0, 1;

	//transformacija koordinata originala i slika u normalizovane koordinate
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

	cout << "nove normalizovane tacke originala:" << endl << endl;
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

	return lastMatrix;
}

int main() {

	//ovde cuvamo tacke na koje kliknemo u interfejsu ili ako smo izabrali drugu opciju tacke koje unesemo rucno
	vector<vector<double>> org_p;

	//broj tacaka koje unosi korisnik, default = 4
	int n = 4;

	//ovde privremeno cuvamo tacke na koje kliknemo u interfejsu 
	vector<double> tmp;

	//visina i sirina slike
	int width = 600;
	int height = 400;
	cout << "Ako zelite da sami unosite originalne tacke unesite broj 1" << endl;
	cout << "Ako pak zelite da birate tacke kao piksele slike koja se ispravlja unesite 2" << endl;
	cout << "Napomena: ako ste izbrali opciju 1 ne dolazi do ispisa slike..." << endl;
	cout << "odgovor: ";
	int odgovor;
	cin >> odgovor;
	
	//load an image

	//kreira se slika koju ucitavamo iz bmp fajla
	//ako zelite da testirate na nekoj svojoj slici - premestite tu sliku u folder gde se nalazi ovaj cpp fajl
	//neka vam slika bude u bmp formatu za svaki slucaj 
	//uvde samo unesite njeno ime
	CImg<unsigned char> image("Perfection.bmp");
	//pomocna slika zato sto na prvoj imam neko obelezavanje izabranih tacaka koje ne zelim da imam na konstruisanoj slici
	CImg<unsigned char> image2("Perfection.bmp");

	//pravimo novu sliku koja je iste visine i sirine kao originalna ona je skroz crna
	CImg<unsigned char> visu(width, height, 1, 3, 0);

	if (odgovor == 2) {

		//svaku sliku koja nam dodje reskaliramo da izgleda interfejs lepse
		image.resize(width, height);
		image2.resize(width, height);

		//pravimo prozore i ucitavamo slike u njih
		CImgDisplay original_picture(image, "Click 4 or more points - press Enter to continue");

		//pomeramo prozore da interfejs izgleda lepse
		original_picture.move(383, 150);


		//dok ne pritisnemo enter ili zatvorimo prozore radi se unos tacaka (samo klikni negde na unetoj slici)
		while (!original_picture.is_keyENTER() && !original_picture.is_closed()) {
			original_picture.wait();
			if (original_picture.button() && original_picture.mouse_y() >= 0) {

				//unos koordinata kliknute tacke
				const int y = original_picture.mouse_y();
				const int x = original_picture.mouse_x();

				//unosimo koordinate kliknutog piksela u privremeni vector
				tmp.push_back(x);
				tmp.push_back(y);
				tmp.push_back(1);
				//click_points.push_back(click_tmp);

				//ovde taj privremeni vektor upisujemo u vector originalnih_tacaka
				org_p.push_back(tmp);
				tmp.clear();


				std::cout << "(" << x << "," << y << ") ="
				<< " R" << (int)image(x, y, 0, 0)
				<< " G" << (int)image(x, y, 0, 1)
				<< " B" << (int)image(x, y, 0, 2) << endl;

				visu(x, y, 0, 0) = (int)image(x, y, 0, 0); //R
				visu(x, y, 0, 1) = (int)image(x, y, 0, 1); //G
				visu(x, y, 0, 2) = (int)image(x, y, 0, 2); //B

				//iscrtavaju se plavi kvadrati tamo gde korisnik klikne na piksel
				for (int i = 0; i < 8; i++)
				{
					for (int j = 0; j < 8; j++) {
						image(x + i, y + j, 0, 0) = 0;   //R
						image(x + i, y + j, 0, 1) = 0;   //G
						image(x + i, y + j, 0, 2) = 255; //B

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
			}
	}

	cout << "originalne tacke:" << endl;
	for (vector<double> p : org_p) {
		for (double e : p) {
			std::cout << e << " ";
		}
		std::cout << endl;
	}
	cout << endl;

	n = org_p.size();

	}
	
	//Unos koordinata ako je izabrana opcija 1
	
	if (odgovor == 1) {
		//ovaj deo ucitava originalne tacke(A, B, C, D) figure koje zelimo da projektujemo

		//pomocni int za ucitavanje
		double x;

		//brojac za pojedinacne tacke, kad prodje 3 broja resetuje se
		int counter = 0;

		//brojac za petlju 
		int brojac = 0;
		
		cout << "unesite broj tacaka" << endl;
		cin >> n;
		cout << "unesite originalne tacke: " << endl;
		while (brojac < (n * 2)) {

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
	}

	//Sve isto ko za originalne tacke - ovde se samo radi o tackama u koje zelimo da projektujemo nase originalne
	// tacke Ap, Bp, Cp, Dp npr.
	vector<vector<double>> dst_p;
	double x2;
	int counter2 = 0, brojac2 = 0;

	std::cout << "Unesite tacke slike: " << endl;

	while (brojac2 < (n*2)) {

		brojac2++;
		std::cin >> x2;
		counter2++;
		tmp.push_back(x2);

		if (counter2 == 2)
		{
			tmp.push_back(1);
			dst_p.push_back(tmp);
			tmp.clear();
			counter2 = 0;
		}
	}
	
	//Naivni algoritam

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

	//iscrtavanje izmenjene slike
	if (odgovor == 2) {

		//drugi prozor gde crtamo projektovanu sliku
		CImgDisplay draw_disp(visu, "The processed image");

		//Prozor sa originalnom slikom da mogu da se vide napravljene razlike
		CImgDisplay original_picture_2(image, "Original picture");

		//pomeramo prozore da interfejs izgleda lepse
		original_picture_2.move(83, 130);
		draw_disp.move(width+83, 130);


		//pomocne promenljive za projektovanje
		MatrixXd TmpOrgPoint;
		TmpOrgPoint.resize(3, 1);

		MatrixXd TmpDstPoint;
		TmpDstPoint.resize(3, 1);

		int cp, rp, cpp, rpp;

		//dva nacina za konstrukciju izmenjene slike - kada iteriramo po pikselima originalne slike i koristimo Q matricu
		//a mozemo i da iteriramo kroz sliku i koristimo matricu Q.inverse()
		//radimo drugi nacin zato sto daje bolje rezultate

		//invertovana Q matrica - matrica dobijena naivnim algoritmom
		Matrix3d Qinv;
		Qinv = Q.inverse();

		for (int r = 0; r < height; r++)
			for (int c = 0; c < width; c++)
			{
				TmpDstPoint << c, r, 1;

				//trzimo koordinate pikela originalne slike koji odgovaraju trenutnom pikselu slike koju konstruisemo
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

				//iscrtavamo piksel na odgovarajucoj lokaciji sa istom bojom kao u originalnoj slici
				visu(c, r, 0, 0) = (int)image2(cp, rp, 0, 0); //R
				visu(c, r, 0, 1) = (int)image2(cp, rp, 0, 1); //G
				visu(c, r, 0, 2) = (int)image2(cp, rp, 0, 2); //B

				visu(c, r, 0, 0) = (int)image2(cpp, rpp, 0, 0); //R
				visu(c, r, 0, 1) = (int)image2(cpp, rpp, 0, 1); //G`	
				visu(c, r, 0, 2) = (int)image2(cpp, rpp, 0, 2); //B

			}
		//prikaz slika u prozorima - dok neko ne ugasi te prozore
		visu.display(draw_disp);
		image.display(original_picture_2);
		while (!draw_disp.is_closed() && !original_picture_2.is_closed()) {
			draw_disp.wait();
		}
	}

	//DLT

		//DLT matrica treba da bude slicna kao Q matrica koju smo dobili iz naivnog algoritma
		Matrix3d V = DLT_ALG(org_p, dst_p, n);
		cout << "Matrica dobijena DLT algoritmom: " << endl << endl;
		cout << V << endl << endl;

		Matrix3d DLT_skalirano;
		//skaliramo dobijenu matricu tako sto je delimo sa V(0,0) i mnozimo sa Q(0, 0)
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
		
		Matrix3d lastMatrix = normalizovani_DLT(org_p, dst_p, n);

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



		//Provera da li idemo iz D u Dp
		MatrixXd Provera4;

		Provera4 = skalirano_normirano_DLT * D;

		cout << "provera da li idemo iz D u Dp sa skaliranom normiranom_DLT matricom:" << endl << endl;
		cout << Provera4 << endl << endl;
		
		//testiranje da li je DLT otporan na transformaciju koordinata

		// da li je matrica P koja se dobije DLT alg ista kao i matrica 
		// Tp.inv * Ptran * T - po svemu bi trebalo, ali...

		//stvaram "nasumicne" transformacije za tacke originala i slike
		Matrix3d testT;

		testT << 0, 1, 2,
				 -1, 0, 3,
				 0, 0, 1;

		Matrix3d testTp;

		testTp << 1, -1, 5,
				 1, 1, -2,
				 0, 0, 1;

		//primena transformacija na nase tacke
		//transformacija koordinata originala i slika u druge koordinate
		
		MatrixXd tmpPoint2;
		MatrixXd tmpNormPoint2;
		tmpPoint2.resize(3, 1);

		vector<vector <double>> ch_org_p;
		vector<double> tmpVec2;
		for (int i = 0; i < n; i++) {
			tmpPoint2 << org_p[i][0], org_p[i][1], org_p[i][2];
			tmpNormPoint2 = testT * tmpPoint2;
			tmpVec2.push_back(tmpNormPoint2(0, 0));
			tmpVec2.push_back(tmpNormPoint2(1, 0));
			tmpVec2.push_back(tmpNormPoint2(2, 0));
			ch_org_p.push_back(tmpVec2);
			tmpVec2.clear();
		}

		cout << "nove transformisane tacke originala:" << endl << endl;
		for (vector<double> aca : ch_org_p) {
			for (double a : aca) {
				cout << a << " ";
			}
			cout << endl;
		}
		cout << endl;

		vector<vector <double>> ch_org_p_s;

		for (int i = 0; i < n; i++) {
			tmpPoint2 << dst_p[i][0], dst_p[i][1], dst_p[i][2];
			tmpNormPoint2 = testTp * tmpPoint2;
			tmpVec2.push_back(tmpNormPoint2(0, 0));
			tmpVec2.push_back(tmpNormPoint2(1, 0));
			tmpVec2.push_back(tmpNormPoint2(2, 0));
			ch_org_p_s.push_back(tmpVec2);
			tmpVec2.clear();
		}
		cout << "nove transformisane tacke slike:" << endl << endl;
		for (vector<double> aca : ch_org_p_s) {
			for (double a : aca) {
				cout << a << " ";
			}
			cout << endl;
		}
		cout << endl;

		//radimo DLT za nove tacke
		Matrix3d DLT_nad_transformisanim_tackama = DLT_ALG(ch_org_p, ch_org_p_s, n);

		cout << "DLT za nove tacke" << endl << endl;
		cout << DLT_nad_transformisanim_tackama << endl << endl;

		//Tp.inv * Ptran * T - racunamo ovo
		Matrix3d actualDLT;

		actualDLT = testTp.inverse() * DLT_nad_transformisanim_tackama * testT;
		cout << "DLT koji bi trebalo da transformise originalne tacke:" << endl << endl;
		cout << actualDLT << endl << endl;

		cout << "Matrica dobijena DLT algoritmom: " << endl << endl;
		cout << V << endl << endl;

		//skaliramo dobijeni algoritam
		Matrix3d skalirano_actualDLT;

		skalirano_actualDLT << actualDLT(0, 0)* Q(0, 0) / actualDLT(0, 0)
			, actualDLT(0, 1)* Q(0, 0) / actualDLT(0, 0)
			, actualDLT(0, 2)* Q(0, 0) / actualDLT(0, 0)
			, actualDLT(1, 0)* Q(0, 0) / actualDLT(0, 0)
			, actualDLT(1, 1)* Q(0, 0) / actualDLT(0, 0)
			, actualDLT(1, 2)* Q(0, 0) / actualDLT(0, 0)
			, actualDLT(2, 0)* Q(0, 0) / actualDLT(0, 0)
			, actualDLT(2, 1)* Q(0, 0) / actualDLT(0, 0)
			, actualDLT(2, 2)* Q(0, 0) / actualDLT(0, 0);

		cout << "Skaliran DLT koji smo dobili iz Tp.inv() * DLT * T " << endl << endl;
		cout << skalirano_actualDLT << endl << endl;

		cout << "Matrica koja se dobije kao rezultat DLT-a na originalnim tackama" << endl << endl;
		cout << DLT_skalirano << endl << endl;
	
		return 0;
}