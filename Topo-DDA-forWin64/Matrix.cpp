#include "definition.h"
//---------------------------------------------------------------------Matrixcd---------------------------------------------------------------------
 Matrixcd::Matrixcd() {
}

 Matrixcd::Matrixcd(int input_size1, int input_size2) {
	m_vec = vector<vector<complex<double>>>(input_size1, vector<complex<double>>(input_size2));
}

 Matrixcd::Matrixcd(const Matrixcd& copy) {
	m_vec = copy.m_vec;
}

 Matrixcd& Matrixcd::operator= (const Matrixcd& input)
{
	// do the copy
	m_vec = input.m_vec;

	// return the existing object so we can chain this operator
	return *this;
}


complex<double>& Matrixcd::operator()(const int pos1, const int pos2) {
	return m_vec[pos1][pos2];
}

//Vector * Scalar

Matrixcd Matrixcd::operator*(const complex<double>& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixcd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] * input;
		}
	}
	return result;
}


Matrixcd Matrixcd::operator*(const double& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixcd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] * input;
		}
	}
	return result;
}


Matrixcd Matrixcd::operator*(const int& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixcd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] * double(input);
		}
	}
	return result;
}

//Vector / Scalar

Matrixcd Matrixcd::operator/(const complex<double>& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixcd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] / input;
		}
	}
	return result;
}


Matrixcd Matrixcd::operator/(const double& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixcd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] / input;
		}
	}
	return result;
}


Matrixcd Matrixcd::operator/(const int& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixcd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] / double(input);
		}
	}
	return result;
}


void Matrixcd::print() {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			cout << m_vec[i][j] << ",";
		}
		cout << endl;
	}
	return;
}


Vectori Matrixcd::getShape() {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Vectori result(2);
	result(0) = length1;
	result(1) = length2;
	return result;
}

//---------------------------------------------------------------------Matrixd---------------------------------------------------------------------
Matrixd::Matrixd() {
}

Matrixd::Matrixd(int input_size1, int input_size2) {
	m_vec = vector<vector<double>>(input_size1, vector<double>(input_size2));
}

Matrixd::Matrixd(const Matrixd& copy) {
	m_vec = copy.m_vec;
}

Matrixd& Matrixd::operator= (const Matrixd& input)
{
	// do the copy
	m_vec = input.m_vec;

	// return the existing object so we can chain this operator
	return *this;
}


double& Matrixd::operator()(const int pos1, const int pos2) {
	return m_vec[pos1][pos2];
}

//Vector * Scalar

Matrixcd Matrixd::operator*(const complex<double>& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixcd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] * input;
		}
	}
	return result;
}


Matrixd Matrixd::operator*(const double& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] * input;
		}
	}
	return result;
}


Matrixd Matrixd::operator*(const int& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] * double(input);
		}
	}
	return result;
}

//Vector / Scalar

Matrixcd Matrixd::operator/(const complex<double>& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixcd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] / input;
		}
	}
	return result;
}


Matrixd Matrixd::operator/(const double& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] / input;
		}
	}
	return result;
}


Matrixd Matrixd::operator/(const int& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] / double(input);
		}
	}
	return result;
}


void Matrixd::print() {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			cout << m_vec[i][j] << ",";
		}
		cout << endl;
	}
	return;
}


Vectori Matrixd::getShape() {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Vectori result(2);
	result(0) = length1;
	result(1) = length2;
	return result;
}

//---------------------------------------------------------------------Matrixi---------------------------------------------------------------------
Matrixi::Matrixi() {
}

Matrixi::Matrixi(int input_size1, int input_size2) {
	m_vec = vector<vector<int>>(input_size1, vector<int>(input_size2));
}

Matrixi::Matrixi(const Matrixi& copy) {
	m_vec = copy.m_vec;
}

Matrixi& Matrixi::operator= (const Matrixi& input)
{
	// do the copy
	m_vec = input.m_vec;

	// return the existing object so we can chain this operator
	return *this;
}


int& Matrixi::operator()(const int pos1, const int pos2) {
	return m_vec[pos1][pos2];
}

//Vector * Scalar

Matrixcd Matrixi::operator*(const complex<double>& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixcd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = double(m_vec[i][j]) * input;
		}
	}
	return result;
}


Matrixd Matrixi::operator*(const double& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = double(m_vec[i][j]) * input;
		}
	}
	return result;
}


Matrixi Matrixi::operator*(const int& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixi result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] * input;
		}
	}
	return result;
}

//Vector / Scalar

Matrixcd Matrixi::operator/(const complex<double>& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixcd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = double(m_vec[i][j]) / input;
		}
	}
	return result;
}


Matrixd Matrixi::operator/(const double& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixd result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = double(m_vec[i][j]) / input;
		}
	}
	return result;
}


Matrixi Matrixi::operator/(const int& input) {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Matrixi result(length1, length2);
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			result(i, j) = m_vec[i][j] / input;
		}
	}
	return result;
}


void Matrixi::print() {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	for (int i = 0; i < length1; i++) {
		for (int j = 0; j < length2; j++) {
			cout << m_vec[i][j] << ",";
		}
		cout << endl;
	}
	return;
}


Vectori Matrixi::getShape() {
	int length1 = m_vec.size();
	int length2 = m_vec[0].size();
	Vectori result(2);
	result(0) = length1;
	result(1) = length2;
	return result;
}

//---------------------------------------------------------------------Overall---------------------------------------------------------------------
Matrixcd matadd(Matrixcd x, Matrixcd y) {
	Vectori shape1 = x.getShape();
	Vectori shape2 = y.getShape();
	if (shape1(0) != shape2(0) || shape1(1) != shape2(1)) {
		cout << "ERROR operator+(Matrixcd& x, Matrixcd& y) : shape not equal" << endl;
		throw 1;
	}
	Matrixcd result(shape1(0), shape1(1));
	for (int i = 0; i < shape1(0); i++) {
		for (int j = 0; j < shape1(1); j++) {
			result(i, j) = x(i, j) + y(i, j);
		}
	}
	return result;
}