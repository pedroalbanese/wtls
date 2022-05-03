// Parameters for the Wireless Transport Layer Security (WAP-WTLS) curves
package wtls

import (
	"crypto/elliptic"
	"errors"
	"math/big"
	"sync"
)

var initonce sync.Once
var p112 *Curve
var p160 *Curve

func initP112() {
	p112 = new(Curve)
	p112.P, _ = new(big.Int).SetString("fffffffffffffffffffffffffde7", 16)
	p112.N, _ = new(big.Int).SetString("0100000000000001ecea551ad837e9", 16)
	p112.B, _ = new(big.Int).SetString("03", 16)
	p112.Gx, _ = new(big.Int).SetString("01", 16)
	p112.Gy, _ = new(big.Int).SetString("02", 16)
	p112.BitSize = 112
}

func P112() elliptic.Curve {
	initonce.Do(initP112)
	return p112
}

func initP160() {
	p160 = new(Curve)
	p160.P, _ = new(big.Int).SetString("fffffffffffffffffffffffffffffffffffc808f", 16)
	p160.N, _ = new(big.Int).SetString("0100000000000000000001cdc98ae0e2de574abf33", 16)
	p160.B, _ = new(big.Int).SetString("03", 16)
	p160.Gx, _ = new(big.Int).SetString("01", 16)
	p160.Gy, _ = new(big.Int).SetString("02", 16)
	p160.BitSize = 160
}

func P160() elliptic.Curve {
	initonce.Do(initP160)
	return p160
}

type Curve struct {
	P       *big.Int // the order of the underlying field
	N       *big.Int // the order of the base point
	B       *big.Int // the constant of the Curve equation
	Gx, Gy  *big.Int // (x,y) of the base point
	BitSize int      // the size of the underlying field
}

func (curve *Curve) Params() *elliptic.CurveParams {
	return &elliptic.CurveParams{
		P:       curve.P,
		N:       curve.N,
		B:       curve.B,
		Gx:      curve.Gx,
		Gy:      curve.Gy,
		BitSize: curve.BitSize,
	}
}

func (curve *Curve) IsOnCurve(x, y *big.Int) bool {
	y2 := new(big.Int).Mul(y, y) 
	y2.Mod(y2, curve.P)          

	x3 := new(big.Int).Mul(x, x) 
	x3.Mul(x3, x)                

	x3.Add(x3, curve.B)
	x3.Mod(x3, curve.P)

	return x3.Cmp(y2) == 0
}

func (curve *Curve) affineFromJacobian(x, y, z *big.Int) (xOut, yOut *big.Int) {
	zinv := new(big.Int).ModInverse(z, curve.P)
	zinvsq := new(big.Int).Mul(zinv, zinv)

	xOut = new(big.Int).Mul(x, zinvsq)
	xOut.Mod(xOut, curve.P)
	zinvsq.Mul(zinvsq, zinv)
	yOut = new(big.Int).Mul(y, zinvsq)
	yOut.Mod(yOut, curve.P)
	return
}

func (curve *Curve) Add(x1, y1, x2, y2 *big.Int) (*big.Int, *big.Int) {
	z := new(big.Int).SetInt64(1)
	return curve.affineFromJacobian(curve.addJacobian(x1, y1, z, x2, y2, z))
}

func (curve *Curve) addJacobian(x1, y1, z1, x2, y2, z2 *big.Int) (*big.Int, *big.Int, *big.Int) {
	z1z1 := new(big.Int).Mul(z1, z1)
	z1z1.Mod(z1z1, curve.P)
	z2z2 := new(big.Int).Mul(z2, z2)
	z2z2.Mod(z2z2, curve.P)

	u1 := new(big.Int).Mul(x1, z2z2)
	u1.Mod(u1, curve.P)
	u2 := new(big.Int).Mul(x2, z1z1)
	u2.Mod(u2, curve.P)
	h := new(big.Int).Sub(u2, u1)
	if h.Sign() == -1 {
		h.Add(h, curve.P)
	}
	i := new(big.Int).Lsh(h, 1)
	i.Mul(i, i)
	j := new(big.Int).Mul(h, i)

	s1 := new(big.Int).Mul(y1, z2)
	s1.Mul(s1, z2z2)
	s1.Mod(s1, curve.P)
	s2 := new(big.Int).Mul(y2, z1)
	s2.Mul(s2, z1z1)
	s2.Mod(s2, curve.P)
	r := new(big.Int).Sub(s2, s1)
	if r.Sign() == -1 {
		r.Add(r, curve.P)
	}
	r.Lsh(r, 1)
	v := new(big.Int).Mul(u1, i)

	x3 := new(big.Int).Set(r)
	x3.Mul(x3, x3)
	x3.Sub(x3, j)
	x3.Sub(x3, v)
	x3.Sub(x3, v)
	x3.Mod(x3, curve.P)

	y3 := new(big.Int).Set(r)
	v.Sub(v, x3)
	y3.Mul(y3, v)
	s1.Mul(s1, j)
	s1.Lsh(s1, 1)
	y3.Sub(y3, s1)
	y3.Mod(y3, curve.P)

	z3 := new(big.Int).Add(z1, z2)
	z3.Mul(z3, z3)
	z3.Sub(z3, z1z1)
	if z3.Sign() == -1 {
		z3.Add(z3, curve.P)
	}
	z3.Sub(z3, z2z2)
	if z3.Sign() == -1 {
		z3.Add(z3, curve.P)
	}
	z3.Mul(z3, h)
	z3.Mod(z3, curve.P)

	return x3, y3, z3
}

func (curve *Curve) Double(x1, y1 *big.Int) (*big.Int, *big.Int) {
	z1 := new(big.Int).SetInt64(1)
	return curve.affineFromJacobian(curve.doubleJacobian(x1, y1, z1))
}

func (curve *Curve) doubleJacobian(x, y, z *big.Int) (*big.Int, *big.Int, *big.Int) {
	a := new(big.Int).Mul(x, x) 
	b := new(big.Int).Mul(y, y) 
	c := new(big.Int).Mul(b, b) 

	d := new(big.Int).Add(x, b) 
	d.Mul(d, d)                 
	d.Sub(d, a)                 
	d.Sub(d, c)                 
	d.Mul(d, big.NewInt(2))     

	e := new(big.Int).Mul(big.NewInt(3), a) 
	f := new(big.Int).Mul(e, e)             

	x3 := new(big.Int).Mul(big.NewInt(2), d) 
	x3.Sub(f, x3)                            
	x3.Mod(x3, curve.P)

	y3 := new(big.Int).Sub(d, x3)                  
	y3.Mul(e, y3)                                  
	y3.Sub(y3, new(big.Int).Mul(big.NewInt(8), c)) 
	y3.Mod(y3, curve.P)

	z3 := new(big.Int).Mul(y, z) 
	z3.Mul(big.NewInt(2), z3)    
	z3.Mod(z3, curve.P)

	return x3, y3, z3
}

func (curve *Curve) ScalarMult(Bx, By *big.Int, k []byte) (*big.Int, *big.Int) {
	Bz := new(big.Int).SetInt64(1)
	x := Bx
	y := By
	z := Bz

	seenFirstTrue := false
	for _, byte := range k {
		for bitNum := 0; bitNum < 8; bitNum++ {
			if seenFirstTrue {
				x, y, z = curve.doubleJacobian(x, y, z)
			}
			if byte&0x80 == 0x80 {
				if !seenFirstTrue {
					seenFirstTrue = true
				} else {
					x, y, z = curve.addJacobian(Bx, By, Bz, x, y, z)
				}
			}
			byte <<= 1
		}
	}

	if !seenFirstTrue {
		return nil, nil
	}

	return curve.affineFromJacobian(x, y, z)
}

func (curve *Curve) ScalarBaseMult(k []byte) (*big.Int, *big.Int) {
	return curve.ScalarMult(curve.Gx, curve.Gy, k)
}

func CompressPoint(curve *Curve, X, Y *big.Int) (cp []byte) {
    return curve.CompressPoint(X, Y)
}

func (curve *Curve) CompressPoint(X, Y *big.Int) (cp []byte) {
	by := new(big.Int).And(Y, big.NewInt(1)).Int64()
	bx := X.Bytes()
	cp = make([]byte, len(bx)+1)
	if by == 1 {
		cp[0] = byte(3)
	} else {
		cp[0] = byte(2)
	}
	copy(cp[1:], bx)

	return
}

func (curve *Curve) DecompressPoint(cp []byte) (X, Y *big.Int, err error) {
	var c int64

	switch cp[0] { 
	case byte(0x03):
		c = 1
		break
	case byte(0x02):
		c = 0
		break
	case byte(0x04):
		X, Y = elliptic.Unmarshal(curve, cp)
		return
	default:
		return nil, nil, errors.New("Not a compressed point. (Invalid Header)")
	}

	byteLen := (curve.Params().BitSize + 7) >> 3
	if len(cp) != 1+byteLen {
		return nil, nil, errors.New("Not a compressed point. (Require 1 + key size)")
	}

	X = new(big.Int).SetBytes(cp[1:])
	Y = new(big.Int)

	Y.Mod(Y.Mul(X, X), curve.P)
	Y.Mod(Y.Mul(Y, X), curve.P)
	Y.Mod(Y.Add(Y, curve.B), curve.P)

	Y = curve.Sqrt(Y)

	if Y.Cmp(big.NewInt(0)) == 0 {
		return nil, nil, errors.New("Not a compressed point. (Not on curve)")
	}

	if c != new(big.Int).And(Y, big.NewInt(1)).Int64() {
		Y.Sub(curve.P, Y)
	}

	return
}

func (curve *Curve) Sqrt(a *big.Int) *big.Int {
	ZERO := big.NewInt(0)
	ONE := big.NewInt(1)
	TWO := big.NewInt(2)
	THREE := big.NewInt(3)
	FOUR := big.NewInt(4)

	p := curve.P
	c := new(big.Int)

	if a.Cmp(ZERO) == 0 {
		return ZERO
	} else if p.Cmp(TWO) == 0 {
		return a.Mod(a,p)
	} else if LegendreSymbol(a, p) != 1 {
		return ZERO
	} else if c.Mod(p, FOUR).Cmp(THREE) == 0 {
		c.Add(p, ONE)
		c.Div(c, FOUR)
		c.Exp(a, c, p)
		return c
	}
	s := new(big.Int)
	s.Sub(p, ONE)

	e := new(big.Int)
	e.Set(ZERO)
	for c.Mod(s, TWO).Cmp(ZERO) == 0 {
		s.Div(s, TWO)
		e.Add(e, ONE)
	}
	n := new(big.Int)
	n.Set(TWO)
	for LegendreSymbol(n, p) != -1 {
		n.Add(n, ONE)
	}
	x := new(big.Int)
	x.Add(s, ONE)
	x.Div(x, TWO)
	x.Exp(a, x, p)
	b := new(big.Int)
	b.Exp(a, s, p)
	g := new(big.Int)
	g.Exp(n, s, p)
	r := new(big.Int)
	r.Set(e)

	t := new(big.Int)
	m := new(big.Int)
	gs := new(big.Int)

	for {
		t.Set(b)
		m.Set(ZERO)

		for ; m.Cmp(r) < 0; m.Add(m, ONE) {
			if t.Cmp(ONE) == 0 {
				break
			}
			t.Exp(t, TWO, p)
		}

		if m.Cmp(ZERO) == 0 {
			return x
		}

		gs.Sub(r, m)
		gs.Sub(gs, ONE)
		gs.Exp(TWO, gs, nil)
		gs.Exp(g, gs, p)

		g.Mod(g.Mul(gs, gs), p)
		x.Mod(x.Mul(x, gs), p)
		b.Mod(b.Mul(b, g), p)
		r.Set(m)
	}
}

func LegendreSymbol(a, p *big.Int) int {
	ZERO := big.NewInt(0)
	ONE := big.NewInt(1)
	TWO := big.NewInt(2)

	ls := new(big.Int).Mod(a, p)

	if ls.Cmp(ZERO) == 0 {
		return 0
	}

	ps := new(big.Int).Sub(p, ONE)

	ls.Div(ps, TWO)
	ls.Exp(a, ls, p)

	if c := ls.Cmp(ps); c == 0 {
		return -1
	}

	return 1
}