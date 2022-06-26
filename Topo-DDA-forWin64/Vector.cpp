#include "definition.h"
//---------------------------------------------------------------------Vectorcd---------------------------------------------------------------------
 Vectorcd::Vectorcd() {
}

 Vectorcd::Vectorcd(int input_size, double initial) {
	m_vec = vector<complex<double>>(input_size, initial);
}
 Vectorcd::Vectorcd(vector<complex<double>> input) {
	 m_vec = input;
 }
 Vectorcd::Vectorcd(const Vectorcd& copy) {
	m_vec = copy.m_vec;
}

 Vectorcd& Vectorcd::operator= (const Vectorcd& input)
{   
	if (this == &input) {
		 return *this;
	}
	// do the copy
	m_vec = input.m_vec;

	// return the existing object so we can chain this operator
	return *this;
}

 Vectorcd& Vectorcd::operator= (Vectorcd&& input)
 {   
	 if (this == &input) {
		 return *this;
	 }
	 // do the move
	 m_vec = move(input.m_vec);

	 // return the existing object so we can chain this operator
	 return *this;
 }

complex<double>& Vectorcd::operator() (const int pos) {
	return m_vec[pos];
}

const complex<double>& Vectorcd::operator() (const int pos) const {
	return m_vec[pos];
}


//Vector / Scalar

Vectorcd Vectorcd::operator/(const complex<double>& input) {
	int length = m_vec.size();
	Vectorcd result(length);
	for (int i = 0; i < length; i++) {
		result(i) = m_vec[i] / input;
	}
	return result;
}


Vectorcd Vectorcd::operator/(const double& input) {
	int length = m_vec.size();
	Vectorcd result(length);
	for (int i = 0; i < length; i++) {
		result(i) = m_vec[i] / input;
	}
	return result;
}


Vectorcd Vectorcd::operator/(const int& input) {
	int length = m_vec.size();
	Vectorcd result(length);
	for (int i = 0; i < length; i++) {
		result(i) = m_vec[i] / double(input);
	}
	return result;
}


//Vector += Vector
Vectorcd& Vectorcd::operator+=(const Vectorcd& rhs) {
	int length = rhs.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectorcd::operator+= : input vec size is different." << endl;
		throw 1;
	}
	for (int i = 0; i < length; i++) {
		m_vec[i] += rhs(i);
	}
	return *this;
}
Vectorcd& Vectorcd::operator+=(const Vectord& rhs) {
	int length = rhs.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectorcd::operator+= : input vec size is different." << endl;
		throw 1;
	}
	for (int i = 0; i < length; i++) {
		m_vec[i] += rhs(i);
	}
	return *this;
}
Vectorcd& Vectorcd::operator+=(const Vectori& rhs) {
	int length = rhs.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectorcd::operator+= : input vec size is different." << endl;
		throw 1;
	}
	for (int i = 0; i < length; i++) {
		m_vec[i] += double(rhs(i));
	}
	return *this;
}
//Vector -= Vector
Vectorcd& Vectorcd::operator-=(const Vectorcd& rhs) {
	int length = rhs.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectorcd::operator-= : input vec size is different." << endl;
		throw 1;
	}
	for (int i = 0; i < length; i++) {
		m_vec[i] -= rhs(i);
	}
	return *this;
}
Vectorcd& Vectorcd::operator-=(const Vectord& rhs) {
	int length = rhs.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectorcd::operator-= : input vec size is different." << endl;
		throw 1;
	}
	for (int i = 0; i < length; i++) {
		m_vec[i] -= rhs(i);
	}
	return *this;
}
Vectorcd& Vectorcd::operator-=(const Vectori& rhs) {
	int length = rhs.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectorcd::operator-= : input vec size is different." << endl;
		throw 1;
	}
	for (int i = 0; i < length; i++) {
		m_vec[i] -= double(rhs(i));
	}
	return *this;
}

//Vector dot* Vector

complex<double> Vectorcd::dot(Vectorcd& input) {
	int length = input.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectorcd::dot : input vec size is different." << endl;
		throw 1;
	}
	complex<double> result = 0.0;
	double a, b, c, d;
	for (int i = 0; i < length; i++) {
		a = m_vec[i].real();
		b = m_vec[i].imag();
		c = input(i).real();
		d = input(i).imag();
		result += (a * c + b * d) + (b * c - a * d) * 1.0i;
	}

	return result;
}


complex<double> Vectorcd::dot(Vectord& input) {
	int length = input.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectorcd::dot : input vec size is different." << endl;
		throw 1;
	}
	complex<double> result = 0.0;
	for (int i = 0; i < length; i++) {
		result += m_vec[i] * input(i);
	}
	return result;
}


complex<double> Vectorcd::dot(Vectori& input) {
	int length = input.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectorcd::dot : input vec size is different." << endl;
		throw 1;
	}
	complex<double> result = 0.0;
	for (int i = 0; i < length; i++) {
		result += m_vec[i] * double(input(i));
	}
	return result;
}


void Vectorcd::print() {
	int length = m_vec.size();
	cout << "(";
	for (int i = 0; i < length; i++) {
		cout << m_vec[i] << ",";
	}
	cout << ")" << endl;
	return;
}


int Vectorcd::size() const {
	int length = m_vec.size();

	return length;
}

double Vectorcd::norm() {
	return sqrt((this->dot(*this)).real());
}


complex<double> Vectorcd::sum() {
	complex<double> result = 0.0;
	for (int i = 0; i < m_vec.size(); i++) {
		result += m_vec[i];
	}
	return result;
}


Vectorcd Vectorcd::vecpow(int input) {
	int length = m_vec.size();
	Vectorcd result(length);
	for (int i = 0; i < length; i++) {
		result(i) = pow(m_vec[i], input);
	}
	return result;
}


Vectorcd Vectorcd::cwiseAbs() {
	int length = m_vec.size();
	Vectorcd result(length);
	for (int i = 0; i < length; i++) {
		result(i) = abs(m_vec[i]);
	}
	return result;
}
//---------------------------------------------------------------------Vectord---------------------------------------------------------------------
Vectord::Vectord() {
}

Vectord::Vectord(int input_size, double initial) {
	m_vec = vector<double>(input_size, initial);
}

Vectord::Vectord(vector<double> input) {
	m_vec = input;
}

Vectord::Vectord(const Vectord& copy) {
	m_vec = copy.m_vec;
}

Vectord& Vectord::operator= (const Vectord& input)
{
	if (this == &input) {
		return *this;
	}
	// do the copy
	m_vec = input.m_vec;

	// return the existing object so we can chain this operator
	return *this;
}

Vectord& Vectord::operator= (Vectord&& input)
{
	if (this == &input) {
		return *this;
	}
	// do the move
	m_vec = move(input.m_vec);

	// return the existing object so we can chain this operator
	return *this;
}

double& Vectord::operator()(const int pos) {
	return m_vec[pos];
}

const double& Vectord::operator()(const int pos) const {
	return m_vec[pos];
}


//Vector / Scalar

Vectorcd Vectord::operator/(const complex<double>& input) {
	int length = m_vec.size();
	Vectorcd result(length);
	for (int i = 0; i < length; i++) {
		result(i) = m_vec[i] / input;
	}
	return result;
}


Vectord Vectord::operator/(const double& input) {
	int length = m_vec.size();
	Vectord result(length);
	for (int i = 0; i < length; i++) {
		result(i) = m_vec[i] / input;
	}
	return result;
}


Vectord Vectord::operator/(const int& input) {
	int length = m_vec.size();
	Vectord result(length);
	for (int i = 0; i < length; i++) {
		result(i) = m_vec[i] / double(input);
	}
	return result;
}


//Vector += Vector
Vectord& Vectord::operator+=(const Vectord& rhs) {
	int length = rhs.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectord::operator+= : input vec size is different." << endl;
		throw 1;
	}
	for (int i = 0; i < length; i++) {
		m_vec[i] += rhs(i);
	}
	return *this;
}
Vectord& Vectord::operator+=(const Vectori& rhs) {
	int length = rhs.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectord::operator+= : input vec size is different." << endl;
		throw 1;
	}
	for (int i = 0; i < length; i++) {
		m_vec[i] += rhs(i);
	}
	return *this;
}
//Vector -= Vector
Vectord& Vectord::operator-=(const Vectord& rhs) {
	int length = rhs.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectord::operator-= : input vec size is different." << endl;
		throw 1;
	}
	for (int i = 0; i < length; i++) {
		m_vec[i] -= rhs(i);
	}
	return *this;
}
Vectord& Vectord::operator-=(const Vectori& rhs) {
	int length = rhs.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectord::operator-= : input vec size is different." << endl;
		throw 1;
	}
	for (int i = 0; i < length; i++) {
		m_vec[i] -= rhs(i);
	}
	return *this;
}


//Vector dot* Vector

complex<double> Vectord::dot(Vectorcd& input) {
	int length = input.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectord::dot : input vec size is different." << endl;
		throw 1;
	}
	complex<double> result = 0.0;
	for (int i = 0; i < length; i++) {
	    result += m_vec[i] * input(i);
	}

	return result;
}


double Vectord::dot(Vectord& input) {
	int length = input.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectord::dot : input vec size is different." << endl;
		throw 1;
	}
	double result = 0.0;
	for (int i = 0; i < length; i++) {
		result += m_vec[i] * input(i);
	}
	return result;
}


double Vectord::dot(Vectori& input) {
	int length = input.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectord::dot : input vec size is different." << endl;
		throw 1;
	}
	double result = 0.0;
	for (int i = 0; i < length; i++) {
		result += m_vec[i] * double(input(i));
	}
	return result;
}


void Vectord::print() {
	int length = m_vec.size();
	cout << "(";
	for (int i = 0; i < length; i++) {
		cout << m_vec[i] << ",";
	}
	cout << ")" << endl;
	return;
}


const int Vectord::size() const {
	int length = m_vec.size();

	return length;
}

double Vectord::norm() {
	return sqrt((this->dot(*this)));
}


double Vectord::sum() {
	double result = 0.0;
	for (int i = 0; i < m_vec.size(); i++) {
		result += m_vec[i];
	}
	return result;
}


Vectord Vectord::vecpow(int input) {
	int length = m_vec.size();
	Vectord result(length);
	for (int i = 0; i < length; i++) {
		result(i) = pow(m_vec[i], input);
	}
	return result;
}


Vectord Vectord::cwiseAbs() {
	int length = m_vec.size();
	Vectord result(length);
	for (int i = 0; i < length; i++) {
		result(i) = abs(m_vec[i]);
	}
	return result;
}

//---------------------------------------------------------------------Vectori---------------------------------------------------------------------
Vectori::Vectori() {
}

Vectori::Vectori(int input_size, double initial) {
	m_vec = vector<int>(input_size, initial);
}

Vectori::Vectori(vector<int> input) {
	m_vec = input;
}

Vectori::Vectori(const Vectori& copy) {
	m_vec = copy.m_vec;
}

Vectori& Vectori::operator= (const Vectori& input)
{
	if (this == &input) {
		return *this;
	}
	// do the copy
	m_vec = input.m_vec;

	// return the existing object so we can chain this operator
	return *this;
}

Vectori& Vectori::operator= (Vectori&& input)
{
	if (this == &input) {
		return *this;
	}
	// do the move
	m_vec = move(input.m_vec);

	// return the existing object so we can chain this operator
	return *this;
}

int& Vectori::operator()(const int pos) {
	return m_vec[pos];
}

const int& Vectori::operator()(const int pos) const {
	return m_vec[pos];
}


//Vector / Scalar

Vectorcd Vectori::operator/(const complex<double>& input) {
	int length = m_vec.size();
	Vectorcd result(length);
	for (int i = 0; i < length; i++) {
		result(i) = double(m_vec[i]) / input;
	}
	return result;
}


Vectord Vectori::operator/(const double& input) {
	int length = m_vec.size();
	Vectord result(length);
	for (int i = 0; i < length; i++) {
		result(i) = double(m_vec[i]) / input;
	}
	return result;
}


Vectori Vectori::operator/(const int& input) {
	int length = m_vec.size();
	Vectori result(length);
	for (int i = 0; i < length; i++) {
		result(i) = m_vec[i] / input;
	}
	return result;
}

//Vector += Vector
Vectori& Vectori::operator+=(const Vectori& rhs) {
	int length = rhs.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectori::operator+= : input vec size is different." << endl;
		throw 1;
	}
	for (int i = 0; i < length; i++) {
		m_vec[i] += rhs(i);
	}
	return *this;
}
//Vector -= Vector
Vectori& Vectori::operator-=(const Vectori& rhs) {
	int length = rhs.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectori::operator-= : input vec size is different." << endl;
		throw 1;
	}
	for (int i = 0; i < length; i++) {
		m_vec[i] -= rhs(i);
	}
	return *this;
}
//Vector dot* Vector

complex<double> Vectori::dot(Vectorcd& input) {
	int length = input.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectori::dot : input vec size is different." << endl;
		throw 1;
	}
	complex<double> result = 0.0;
	for (int i = 0; i < length; i++) {
		result += double(m_vec[i]) * input(i);
	}

	return result;
}


double Vectori::dot(Vectord& input) {
	int length = input.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectori::dot : input vec size is different." << endl;
		throw 1;
	}
	double result = 0.0;
	for (int i = 0; i < length; i++) {
		result += double(m_vec[i]) * input(i);
	}
	return result;
}


int Vectori::dot(Vectori& input) {
	int length = input.size();
	if (length != m_vec.size()) {
		cout << "ERROR Vectori::dot : input vec size is different." << endl;
		throw 1;
	}
	int result = 0;
	for (int i = 0; i < length; i++) {
		result += m_vec[i] * input(i);
	}
	return result;
}


void Vectori::print() {
	int length = m_vec.size();
	cout << "(";
	for (int i = 0; i < length; i++) {
		cout << m_vec[i] << ",";
	}
	cout << ")" << endl;
	return;
}


int Vectori::size() const {
	int length = m_vec.size();

	return length;
}

double Vectori::norm() {
	return sqrt((this->dot(*this)));
}


int Vectori::sum() {
	int result = 0;
	for (int i = 0; i < m_vec.size(); i++) {
		result += m_vec[i];
	}
	return result;
}


Vectori Vectori::vecpow(int input) {
	int length = m_vec.size();
	Vectori result(length);
	for (int i = 0; i < length; i++) {
		result(i) = pow(m_vec[i], input);
	}
	return result;
}


Vectori Vectori::cwiseAbs() {
	int length = m_vec.size();
	Vectori result(length);
	for (int i = 0; i < length; i++) {
		result(i) = abs(m_vec[i]);
	}
	return result;
}

//---------------------------------------------------------------------Overall---------------------------------------------------------------------
Vectorcd operator+(const Vectorcd &x, const Vectorcd &y) {
	int length1 = x.size();
	int length2 = y.size();
	if (length1 != length2) {
		cout << "ERROR operator+ : two vec length different" << endl;
	}
	Vectorcd result(length1);
	for (int i = 0; i < length1; i++) {
		result(i) = x(i) + y(i);
	}
	return result;
}
Vectorcd operator+(const Vectord& x, const Vectorcd& y) {
	int length1 = x.size();
	int length2 = y.size();
	if (length1 != length2) {
		cout << "ERROR operator+ : two vec length different" << endl;
	}
	Vectorcd result(length1);
	for (int i = 0; i < length1; i++) {
		result(i) = x(i) + y(i);
	}
	return result;
}
Vectorcd operator+(const Vectorcd& x, const Vectord& y) {
	return y + x;
}
Vectord operator+(const Vectord& x, const Vectord& y) {
	int length1 = x.size();
	int length2 = y.size();
	if (length1 != length2) {
		cout << "ERROR operator+ : two vec length different" << endl;
	}
	Vectord result(length1);
	for (int i = 0; i < length1; i++) {
		result(i) = x(i) + y(i);
	}
	return result;
}




Vectorcd operator-(const Vectorcd& x, const Vectorcd& y) {
	int length1 = x.size();
	int length2 = y.size();
	if (length1 != length2) {
		cout << "ERROR operator- : two vec length different" << endl;
	}
	Vectorcd result(length1);
	for (int i = 0; i < length1; i++) {
		result(i) = x(i) - y(i);
	}
	return result;
}
Vectorcd operator-(const Vectord& x, const Vectorcd& y) {
	int length1 = x.size();
	int length2 = y.size();
	if (length1 != length2) {
		cout << "ERROR operator- : two vec length different" << endl;
	}
	Vectorcd result(length1);
	for (int i = 0; i < length1; i++) {
		result(i) = x(i) - y(i);
	}
	return result;
}
Vectorcd operator-(const Vectorcd& x, const Vectord& y) {
	int length1 = x.size();
	int length2 = y.size();
	if (length1 != length2) {
		cout << "ERROR operator- : two vec length different" << endl;
	}
	Vectorcd result(length1);
	for (int i = 0; i < length1; i++) {
		result(i) = x(i) - y(i);
	}
	return result;
}
Vectord operator-(const Vectord& x, const Vectord& y) {
	int length1 = x.size();
	int length2 = y.size();
	if (length1 != length2) {
		cout << "ERROR operator- : two vec length different" << endl;
	}
	Vectord result(length1);
	for (int i = 0; i < length1; i++) {
		result(i) = x(i) - y(i);
	}
	return result;
}

complex<double> vecmean(const Vectorcd& x) {
	int length = x.size();
	complex<double> result = 0.0;
	for (int i = 0; i < length; i++) {
		result += x(i);
	}
	result = result / double(length);
	return result;
}

double vecmean(const Vectord& x) {
	int length = x.size();
	double result = 0.0;
	for (int i = 0; i < length; i++) {
		result += x(i);
	}
	result = result / length;
	return result;
}

double vecmean(const Vectori& x) {
	int length = x.size();
	double result = 0.0;
	for (int i = 0; i < length; i++) {
		result += x(i);
	}
	result = result / length;
	return result;
}

void vectofile(ofstream& fout, Vectorcd object) {
	for (int i = 0; i < object.size(); i++) {
		fout << object(i) << endl;
	}
	return;
}

void vectofile(ofstream& fout, Vectord object) {
	for (int i = 0; i < object.size(); i++) {
		fout << object(i) << endl;
	}
	return;
}

void vectofile(ofstream& fout, Vectori object) {
	for (int i = 0; i < object.size(); i++) {
		fout << object(i) << endl;
	}
	return;
}

//Vector * Scalar

Vectorcd operator*(const Vectorcd& x, const complex<double>& input) {
	int length = x.size();
	Vectorcd result(length);
	for (int i = 0; i < length; i++) {
		result(i) = x(i) * input;
	}
	return result;
}

Vectorcd operator*(const complex<double>& input, const Vectorcd& x) {
	return x * input;
}


Vectorcd operator*(const Vectorcd& x, const double& input) {
	int length = x.size();
	Vectorcd result(length);
	for (int i = 0; i < length; i++) {
		result(i) = x(i) * input;
	}
	return result;
}

Vectorcd operator*(const double& input, const Vectorcd& x) {
	return x * input;
}

Vectorcd operator*(const Vectorcd& x, const int& input) {
	int length = x.size();
	Vectorcd result(length);
	for (int i = 0; i < length; i++) {
		result(i) = x(i) * double(input);
	}
	return result;
}

Vectorcd operator*(const int& input, const Vectorcd& x) {
	return x * input;
}

//Vector * Scalar

Vectorcd operator*(const Vectord& x, const complex<double>& input) {
	int length = x.size();
	Vectorcd result(length);
	for (int i = 0; i < length; i++) {
		result(i) = x(i) * input;
	}
	return result;
}

Vectorcd operator*(const complex<double>& input, const Vectord& x) {
	return x * input;
}

Vectord operator*(const Vectord& x, const double& input) {
	int length = x.size();
	Vectord result(length);
	for (int i = 0; i < length; i++) {
		result(i) = x(i) * input;
	}
	return result;
}

Vectord operator*(const double& input, const Vectord& x) {
	return x * input;
}

Vectord operator*(const Vectord& x, const int& input) {
	int length = x.size();
	Vectord result(length);
	for (int i = 0; i < length; i++) {
		result(i) = x(i) * double(input);
	}
	return result;
}

Vectord operator*(const int& input, const Vectord& x) {
	return x * input;
}

//Vector * Scalar

Vectorcd operator*(const Vectori& x, const complex<double>& input) {
	int length = x.size();
	Vectorcd result(length);
	for (int i = 0; i < length; i++) {
		result(i) = double(x(i)) * input;
	}
	return result;
}

Vectorcd operator*(const complex<double>& input, const Vectori& x) {
	return x * input;
}

Vectord operator*(const Vectori& x, const double& input) {
	int length = x.size();
	Vectord result(length);
	for (int i = 0; i < length; i++) {
		result(i) = double(x(i)) * input;
	}
	return result;
}

Vectord operator*(const double& input, const Vectori& x) {
	return x * input;
}

Vectori operator*(const Vectori& x, const int& input) {
	int length = x.size();
	Vectori result(length);
	for (int i = 0; i < length; i++) {
		result(i) = x(i) * input;
	}
	return result;
}

Vectori operator*(const int& input, const Vectori& x) {
	return x * input;
}