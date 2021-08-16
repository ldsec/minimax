package main

import(
	"math/big"
	"fmt"
	"math"
)

var pi = "3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679821480865132823066470938446095505822317253594081284811174502841027019385211055596446229489549303819644288109756659334461284756482337867831652712019091456485669234603486104543266482133936072602491412737245870066063155881748815209209628292540917153643678925903600113305305488204665213841469519415116094330572703657595919530921861173819326117931051185480744623799627495673518857527248912279381830119491298336733624406566430860213949463952247371907021798609437027705392171762931767523846748184676694051320005681271452635608277857713427577896091736371787214684409012249534301465495853710507922796892589235420199561121290219608640344181598136297747713099605187072113499999983729780499510597317328160963185950244594553469083026425223082533446850352619311881710100031378387528865875332083814206171776691473035982534904287554687311595628638823537875937519577818577805321712268066130019278766111959092164201989"
var precision uint = 128

type Function interface{
	Eval(x *big.Float) (y *big.Float)
}

type ErrorFunction struct{
	f Function
	a, b *big.Float
	poly []*big.Float
}

func (f *ErrorFunction) Eval(x *big.Float) (y *big.Float){
	y = ChebyshevEval(f.poly, x, f.a, f.b)
	//fmt.Println(x, y)
	y.Sub(y,f.f.Eval(x))
	return
}

type TrigoFunc struct{
}

func (t *TrigoFunc) Eval(x *big.Float) (y *big.Float){
	z := new(big.Float).Set(x)
	z.Neg(z)
	z = Exp(z, precision)
	z.Add(z, NewFloat(1))
	y = NewFloat(1)
	y.Quo(y, z)
	return
}


// NewFloat creates a new big.Float element with 1000 bits of precision
func NewFloat(x float64) (y *big.Float) {
	y = new(big.Float)
	y.SetPrec(precision) // log2 precision
	y.SetFloat64(x)
	return
}

func Logarithm(y *big.Float, precision uint) (ln *big.Float){
	x := new(big.Float).Set(y)
	tmp := new(big.Float).Set(y)
	x.Sub(x, NewFloat(1))
	tmp.Add(tmp, NewFloat(1))
	x.Quo(x, tmp)
	z := new(big.Float).Mul(x, x)
	ln = NewFloat(0)

	for k := uint(1); k < precision; k+=2{
		tmp.Add(x, x)
		tmp.Quo(tmp, NewFloat(float64(k)))
		ln.Add(ln, tmp)
		x.Mul(x, z)
	} 
	return
}

func Exp(x *big.Float, precision uint) (exp *big.Float){
	two := NewFloat(2)

	exp = new(big.Float).Set(x)
	exp.SetPrec(2*precision) // log2 precision
	for i := uint(0); i < precision; i++{
		exp.Quo(exp, two)
	}

	exp.Add(exp, NewFloat(1))

	for i := uint(0); i < precision; i++{
		exp.Mul(exp, exp)
	}

	exp.SetPrec(precision)

	return
}

func SinH(x *big.Float, precision uint) (sinh *big.Float){
	sinh = new(big.Float).Set(x)
	sinh.Add(sinh, sinh)
	sinh.Neg(sinh)
	sinh = Exp(sinh, precision)
	sinh.Neg(sinh)
	sinh.Add(sinh, NewFloat(1))
	tmp := new(big.Float).Set(x)
	tmp.Neg(tmp)
	tmp = Exp(tmp, precision)
	tmp.Add(tmp, tmp)
	sinh.Quo(sinh, tmp)
	return
}

func TanH(x *big.Float, precision uint) (tanh *big.Float){
	tanh = new(big.Float).Set(x)
	tanh.Add(tanh, tanh)
	tanh = Exp(tanh, precision)
	tmp := new(big.Float).Set(tanh)
	tmp.Add(tmp, NewFloat(1))
	tanh.Sub(tanh, NewFloat(1))
	tanh.Quo(tanh, tmp)
	return
}

func ArithmeticMean(x, y *big.Float, k uint) (*big.Float){
	a := new(big.Float).Set(x)
	g := new(big.Float).Set(y)
	tmp := new(big.Float)
	half := NewFloat(0.5)

	for i := 0; i < int(math.Log2(float64(precision))); i++{
		tmp.Mul(a, g)
		a.Add(a, g)
		a.Mul(a, half)
		g.Sqrt(tmp)
	}

	return a
}

// Cos is an iterative arbitrary precision computation of Cos(x)
// Iterative process with an error of ~10^{−0.60206*k} after k iterations.
// ref : Johansson, B. Tomas, An elementary algorithm to evaluate trigonometric functions to high precision, 2018
func Cos(x *big.Float, k uint) (cosx *big.Float) {
	tmp := new(big.Float)

	t := NewFloat(0.5)
	half := new(big.Float).Copy(t)

	for i := uint(1); i < k-1; i++ {
		t.Mul(t, half)
	}

	s := new(big.Float).Mul(x, t)
	s.Mul(s, x)
	s.Mul(s, t)

	four := NewFloat(4.0)

	for i := uint(1); i < k; i++ {
		tmp.Sub(four, s)
		s.Mul(s, tmp)
	}

	cosx = new(big.Float).Quo(s, NewFloat(2.0))
	cosx.Sub(NewFloat(1.0), cosx)
	return

}

func ChebyshevNodes(n int, a, b *big.Float) (nodes []*big.Float) {

	var PiOverN = new(big.Float)
	PiOverN.SetPrec(precision)
	PiOverN.SetString(pi)
	PiOverN.Quo(PiOverN, NewFloat(float64(n-1)))

	nodes = make([]*big.Float, n)

	x := new(big.Float).Add(b, a)
	y := new(big.Float).Sub(b, a)

	two := NewFloat(2)

	x.Quo(x, two)
	y.Quo(y, two)

	for i := 0; i < n; i++{
		nodes[i] = NewFloat(float64(n-i-1))
		nodes[i].Mul(nodes[i], PiOverN)
		nodes[i] = Cos(nodes[i], precision)
		nodes[i].Mul(nodes[i], y)
		nodes[i].Add(nodes[i], x)
	}

	return
}

func ChebyshevBasisInPlace(deg int, x, a, b *big.Float, poly []*big.Float) {
	two := NewFloat(2)

	var tmp, u = new(big.Float), new(big.Float)
	var T, Tprev, Tnext = new(big.Float), new(big.Float), new(big.Float)

	// u = (2*x - (a+b))/(b-a)
	u.Set(x)
	u.Mul(u, two)
	u.Sub(u, a)
	u.Sub(u, b)
	tmp.Set(b)
	tmp.Sub(tmp, a)
	u.Quo(u, tmp)

	Tprev.SetPrec(precision)
	Tprev.SetFloat64(1)
	T.Set(u)
	poly[0].Set(Tprev)

	for i := 1; i < deg; i++{
		Tnext.Mul(two, u)
		Tnext.Mul(Tnext, T)
		Tnext.Sub(Tnext, Tprev)
		Tprev.Set(T)
		T.Set(Tnext)
		poly[i].Set(Tprev)
	}
}

func ChebyshevEval(poly []*big.Float, x, a, b *big.Float) (y *big.Float){

	two := NewFloat(2)
	var tmp, u = new(big.Float), new(big.Float)
	var T, Tprev, Tnext = new(big.Float), new(big.Float), new(big.Float)

	// u = (2*x - (a+b))/(b-a)
	u.Set(x)
	u.Mul(u, two)
	u.Sub(u, a)
	u.Sub(u, b)
	tmp.Set(b)
	tmp.Sub(tmp, a)
	u.Quo(u, tmp)

	Tprev.SetPrec(precision)
	Tprev.SetFloat64(1)
	T.Set(u)
	y = new(big.Float).Set(poly[0])

	for i := 1; i < len(poly); i++{
		y.Add(y, tmp.Mul(T, poly[i]))
		Tnext.Mul(two, u)
		Tnext.Mul(Tnext, T)
		Tnext.Sub(Tnext, Tprev)
		Tprev.Set(T)
		T.Set(Tnext)
	}

	return 
}

func ChebyshevApproximate(f Function, deg int, a, b *big.Float) (coeffs []*big.Float){

	nodes := ChebyshevNodes(deg+1, a, b)

	coeffs = make([]*big.Float, len(nodes))
	for i := range nodes{
		coeffs[i] = new(big.Float)
	}

	two := NewFloat(2)

	var tmp, u = new(big.Float), new(big.Float)
	var T, Tprev, Tnext = new(big.Float), new(big.Float), new(big.Float)
	var y *big.Float

	for i := range nodes{

		y = f.Eval(nodes[i])

		u.Set(nodes[i])
		u.Mul(u, two)
		u.Sub(u, a)
		u.Sub(u, b)
		tmp.Set(b)
		tmp.Sub(tmp, a)
		u.Quo(u, tmp)

		Tprev = NewFloat(1)
		T.Set(u)

		for j := range nodes{
			coeffs[j].Add(coeffs[j], tmp.Mul(y, Tprev))
			Tnext.Mul(two, u)
			Tnext.Mul(Tnext, T)
			Tnext.Sub(Tnext, Tprev)
			Tprev.Set(T)
			T.Set(Tnext)
		}
	}

	tmp = NewFloat(float64(len(nodes)))
	coeffs[0].Quo(coeffs[0], tmp)
	tmp.Quo(tmp, two)
	for i := 1; i < len(nodes); i++{
		coeffs[i].Quo(coeffs[i], tmp)
	}

	return
}

// FindRoot finds a root in the interval a, b
// Expect the function to be strictly increasing or decreasing
// in the interval a, b (hence have only a single root)
func FindRoot(f Function, a, b *big.Float) (root *big.Float){
	two := NewFloat(2)
	atmp := new(big.Float).Set(a)
	btmp := new(big.Float).Set(b)
	root = new(big.Float)

	signRight := f.Eval(b).Sign()

	if f.Eval(a).Sign() == signRight{
		panic("Ill conditionned error function")
	}

	for i := uint(0); i < precision; i++{
		root.Add(atmp, btmp)
		root.Quo(root, two)

		if f.Eval(root).Sign() == signRight{
			btmp.Set(root)
		}else{
			atmp.Set(root)
		}
	}

	return root
}

// FindMaxErr finds a maximum between two roots of a function
func FindMaxErr(f Function, a, b *big.Float) (x, maxErr *big.Float){
	two := NewFloat(2)
	atmp := new(big.Float).Set(a)
	btmp := new(big.Float).Set(b)
	tmp :=  new(big.Float)
	x = new(big.Float)
	for i := uint(0); i < precision; i++{
		x.Add(atmp, btmp)
		x.Quo(x, two)

		err0 := f.Eval(x)
		err0.Abs(err0)

		tmp.Add(atmp, x)
		tmp.Quo(tmp, two)

		err1 := f.Eval(tmp)
		err1.Abs(err1)

		if err0.Cmp(err1) == -1{
			btmp.Set(x)
		}else{
			atmp.Set(x)
		}
	}

	return x, f.Eval(x)
}

// SolveLinearSystemInPlace solves the linear system with Gaussian elimination
func SolveLinearSystemInPlace(matrix [][]*big.Float, vector []*big.Float){

	n, m := len(matrix), len(matrix[0])

	var tmp = new(big.Float)
	for i := 0; i < n; i++{

		a := matrix[i][i]

		vector[i].Quo(vector[i], a)

		for j := m-1; j >= i; j--{
			b := matrix[i][j]
			b.Quo(b, a)
		}

		for j := i+1; j < m; j++{
			c := matrix[j][i]
			vector[j].Sub(vector[j], tmp.Mul(vector[i], c))
			for k := m-1; k >= i; k--{
				matrix[j][k].Sub(matrix[j][k], tmp.Mul(matrix[i][k], c))
			}
		}
	}

	for i := m-1; i > 0; i--{
		c := vector[i]
		for j := i-1; j >=0; j--{
			vector[j].Sub(vector[j], tmp.Mul(matrix[j][i], c))
		}
	}
}

func MinimaxPolynomial(f Function, deg int, a, b *big.Float, tolerance float64) (poly []float64){

	delta := NewFloat(tolerance)

	var Err *big.Float
	minErr := new(big.Float)
	maxErr := new(big.Float)
	normErr := new(big.Float)
	
	matrix := make([][]*big.Float, deg+2)
	vector := make([]*big.Float, deg+2)
	roots := make([]*big.Float, deg+1)

	for i := 0; i < deg+2; i++{
		matrix[i] = make([]*big.Float, deg+2)
		for j := 0; j < deg+1; j++{
			matrix[i][j] = new(big.Float)
		}
	}

	nodes := ChebyshevNodes(deg+2, a, b)

	for {

		for j := 0; j < deg+2; j++{
			ChebyshevBasisInPlace(deg+1, nodes[j], a, b, matrix[j])

			if j&1 == 0{
				matrix[j][deg+1] = NewFloat(-1)
			}else{
				matrix[j][deg+1] = NewFloat(1)
			}

			vector[j] = f.Eval(nodes[j])
		}

		SolveLinearSystemInPlace(matrix, vector)

		fErr := &ErrorFunction{f, a, b, vector[:deg+1]}

		for j := 0; j < deg+1; j++{
			roots[j] = FindRoot(fErr, nodes[j], nodes[j+1])
		}

		var left, right *big.Float

		
		nodes[0], Err = FindMaxErr(fErr, a, roots[0])
		Err.Abs(Err)
		minErr.Set(Err)
		maxErr.Set(Err)
		for j := 1; j < deg+1; j++{

			if j == deg+1{
				left, right = roots[j-1], b
			}else{
				left, right = roots[j-1], roots[j]
			}

			nodes[j], Err = FindMaxErr(fErr, left, right)

			Err.Abs(Err)

			if j > 0{

				if minErr.Cmp(Err) == 1{
					minErr.Set(Err)
				}

				if maxErr.Cmp(Err) == -1{
					maxErr.Set(Err)
				}
			}
		}

		normErr.Sub(maxErr, minErr)
		normErr.Quo(normErr, minErr)

		fmt.Println("diff err :", normErr)

		if normErr.Cmp(delta) < 1{
			poly = make([]float64, deg+1)
			for i := range poly{
				poly[i], _ = vector[i].Float64()
			}
			return
		} 
	}
}



func main(){

	a := NewFloat(-25)
	b := NewFloat(25)
	deg := 32
	
	poly := MinimaxPolynomial(new(TrigoFunc), deg, a, b, 1e-15)
	fmt.Println()

	fmt.Printf("[")
	for i := range poly{
		fmt.Printf("%0.15f,\n", poly[i])
	}
	fmt.Printf("]\n")
}