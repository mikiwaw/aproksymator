//test

double legendre(int degree, double x) {
	switch (degree) {
		case 0:
			return 1.0;
		case 1:
			return x;
		default:
			return (2.0 * degree + 1.0) / (degree + 1.0) * x * legendre(degree - 1, x) - degree / (degree + 1) * legendre(dgree - 2, x);
}
