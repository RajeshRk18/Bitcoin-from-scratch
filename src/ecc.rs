use crate::fieldelement::*;
use crate::traits::{IsEllipticCurve, Secp256k1Curve};
use crate::utils::*;
use lazy_static::lazy_static;
use num_bigint::BigUint;

lazy_static! {
    pub static ref ORDER: FieldElement = FieldElement::new(
        BigUint::parse_bytes(
            "115792089237316195423570985008687907852837564279074904382605163141518161494337"
                .as_bytes(),
            10
        )
        .unwrap()
    );
    pub static ref PRIME: FieldElement = FieldElement::new(
        BigUint::parse_bytes(
            "115792089237316195423570985008687907853269984665640564039457584007908834671663"
                .as_bytes(),
            10
        )
        .unwrap()
    );
    pub static ref G_X: FieldElement = FieldElement::new(
        BigUint::parse_bytes(
            "55066263022277343669578718895168534326250603453777594175500187360389116729240"
                .as_bytes(),
            10
        )
        .unwrap()
    );
    pub static ref G_Y: FieldElement = FieldElement::new(
        BigUint::parse_bytes(
            "32670510020758816978083085130507043184471273380659243275938904335757337482424"
                .as_bytes(),
            10
        )
        .unwrap()
    );

    pub static ref PRECOMPUTED_TABLE: Vec<(FieldElement, FieldElement, FieldElement)> = get_precompute();
}

pub struct Secp256k1;

impl Secp256k1 {

    pub fn mul(&self, scalar: BigUint) -> <Self as IsEllipticCurve>::Affine {
        let proj_result = self.scalar_gen_mul(scalar);
        let proj_result_to_affine: (FieldElement, FieldElement) = self.to_affine(proj_result);

        proj_result_to_affine
    }

    pub fn get_precomputed_table() -> Vec<<Self as IsEllipticCurve>::Jacobian> {
        PRECOMPUTED_TABLE.to_vec()
    }

}

impl Secp256k1Curve for Secp256k1 {
    fn a() -> FieldElement {
        felt_from_uint(0)
    }

    fn b() -> FieldElement {
        felt_from_uint(7)
    }

    fn base_field_order() -> FieldElement {
        ORDER.clone()
    }

    fn scalar_field_order() -> FieldElement {
        PRIME.clone()
    }

    fn generator() -> Self::Affine {
        (G_X.clone(), G_Y.clone())
    }

    fn identity() -> Self::Jacobian {
        (Self::zero(), Self::zero(), Self::one())
    }

    fn one() -> FieldElement {
        felt_from_uint(1)
    }

    fn zero() -> FieldElement {
        felt_from_uint(0)
    }
}

impl IsEllipticCurve for Secp256k1 {
    type Affine = (FieldElement, FieldElement);
    type Jacobian = (FieldElement, FieldElement, FieldElement);
    type Scalar = BigUint;

    #[inline(always)]
    fn is_valid(&self, point: &Self::Affine) -> bool {
        let lhs = &point.1 * &point.1;
        let rhs: FieldElement = &point.0 * &point.0 * &point.0 + felt_from_uint(7);
        if lhs != rhs {
            return false;
        }
        true
    }

    fn is_indiscriminant() -> bool {
        let op = -felt_from_uint(16) * (felt_from_uint(4) * Self::a() * Self::a() * Self::a()
            + felt_from_uint(27) * Self::b() * Self::b());
        if op != Self::zero() {
            return true;
        }
        false
    }

    fn is_negative(point1: &Self::Jacobian, point2: &Self::Jacobian) -> bool {
        if (&point1.0 == &point2.0)
            && (&point1.2 == &point2.2)
            && (point1.1.clone() == -point2.1.clone())
        {
            return true;
        }
        false
    }

    fn is_identity(p: &Self::Jacobian) -> bool {
        let mask = &p.0 == &0 && &p.1 == &0 && &p.2 == &1;
        mask
    }

    #[inline(always)]
    fn to_jacobian(&self, point: Self::Affine) -> Self::Jacobian {
        (point.0, point.1, Self::one())
    }

    #[inline(always)]
    fn to_affine(&self, point: Self::Jacobian) -> Self::Affine {
        let x = point.0.clone();
        let y = point.1.clone();
        let z = point.2.clone();

        if z == Self::zero() {
            return (Self::zero(), Self::zero());
        } else {
            let z_inv =  Self::one() / z;
            let z_inv_squared = &z_inv * &z_inv;
            let z_inv_cubed = &z_inv_squared * &z_inv;

            let x = x * z_inv_squared;
            let y = y * z_inv_cubed;

            return (x, y);
        }
    }

    #[inline(always)]
    fn scalar_gen_mul(&self, scalar: BigUint) -> Self::Jacobian {
        assert!(
            scalar != biguint_from_uint(0 as u32) && FieldElement::from(scalar.clone()) < ORDER.clone()
        );
        let mut scalar = scalar;
        let precomputes = Self::get_precomputed_table();

        let mut result: (FieldElement, FieldElement, FieldElement) = Self::identity();
        let mut dummy_point: (FieldElement, FieldElement, FieldElement) = Self::identity();

        for i in 0..256 {
            if &scalar & biguint_from_uint(1 as u32) == biguint_from_uint(1 as u32) {
                result = Self::point_add(&result, &precomputes[i]);
                dbg!(&result);
            } else {
                // Dummy op to make it constant time.
                dummy_point = Self::point_add(&dummy_point, &precomputes[i]);
            }
            scalar >>= 1;
        }

        result
    }

    #[inline(always)]
    fn point_add(point1: &Self::Jacobian, point2: &Self::Jacobian) -> Self::Jacobian {
        // https://www.hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/addition/add-2007-bl.op3
        /*
        Z1Z1 = Z1^2
        Z2Z2 = Z2^2
        U1 = X1*Z2Z2
        U2 = X2*Z1Z1
        t0 = Z2*Z2Z2
        S1 = Y1*t0
        t1 = Z1*Z1Z1
        S2 = Y2*t1
        H = U2-U1
        t2 = 2*H
        I = t2^2
        J = H*I
        t3 = S2-S1
        r = 2*t3
        V = U1*I
        t4 = r^2
        t5 = 2*V
        t6 = t4-J
        X3 = t6-t5
        t7 = V-X3
        t8 = S1*J
        t9 = 2*t8
        t10 = r*t7
        Y3 = t10-t9
        t11 = Z1+Z2
        t12 = t11^2
        t13 = t12-Z1Z1
        t14 = t13-Z2Z2
        Z3 = t14*H
        */
        let (x1, y1, z1) = point1.clone();
        let (x2, y2, z2) = point2.clone();

        let z1_squared = &z1 * &z1;
        let z2_squared = &z2 * &z2;
        let t1 = &z1_squared * &z1;
        let t0 = &z2_squared * &z2;
        let u1 = &x1 * &z2_squared;
        let u2 = &x2 * &z1_squared;

        let s1 = &t0 * &y1;
        let s2 = &t1 * &y2;
        let h = &u2 - &u1;
        let t2 = felt_from_uint(2) * &h;
        let i = &t2 * &t2;
        let j = &h * &i;
        let t3 = &s2 - &s1;
        let r = felt_from_uint(2) * &t3;
        let v = &u1 * &i;
        let t4 = &r * &r;
        let t5 = felt_from_uint(2) * &v;
        let t6 = &t4 - &j;
        let x3 = &t6 - &t5;
        let t7 = &v - &x3;
        let t8 = &s1 * &j;
        let t9 = felt_from_uint(2) * &t8;
        let t10 = &r * &t7;
        let y3 = &t10 - &t9;
        let t11 = &z1 + &z2;
        let t12 = &t11 * &t11;
        let t13 = &t12 * &z1_squared;
        let t14 = &t13 * &z2_squared;
        let z3 = t14 * &h;

        /*let point1_is_zero = x1 == 0 || y1 == 0;
        let point2_is_zero = x2 == 0 || y2 == 0;*/
        let is_negative = Self::is_negative(&point1, &point2);
        let point1_is_infinity = Self::is_identity(&point1);
        let point2_is_infinity = Self::is_identity(&point2);
        let h_is_zero = h == 0;
        let t3_is_zero = t3 == 0;

        if point1_is_infinity {
            // Dummy op to make it constant time.
            #[allow(unused_variables)]
            let dummy_op = Self::point_double(&point1);
            return point2.clone();
        }

        if point2_is_infinity {
            #[allow(unused_variables)]
            let dummy_op = Self::point_double(&point1);
            return point1.clone();
        }

        if is_negative {
            return Self::point_double(&point1);
        }

        if h_is_zero {
            if t3_is_zero {
                return Self::point_double(&point1);
            } else {
                #[allow(unused_variables)]
                let dummy_op = Self::point_double(&point1);
                return Self::identity();
            }
        } else {
            #[allow(unused_variables)]
            let dummy_op = Self::point_double(&point1);
            return (x3, y3, z3);
        }
    }

    #[inline(always)]
    fn point_double(point: &Self::Jacobian) -> Self::Jacobian {
        // https://www.hyperelliptic.org/EFD/g1p/auto-code/shortw/jacobian-0/doubling/dbl-2007-bl.op3

        /*
        XX = X1^2
        YY = Y1^2
        YYYY = YY^2
        ZZ = Z1^2
        t0 = X1+YY
        t1 = t0^2
        t2 = t1-XX
        t3 = t2-YYYY
        S = 2*t3
        t4 = ZZ^2
        t5 = a*t4
        t6 = 3*XX
        M = t6+t5
        t7 = M^2
        t8 = 2*S
        T = t7-t8
        X3 = T
        t9 = S-T
        t10 = 8*YYYY
        t11 = M*t9
        Y3 = t11-t10
        t12 = Y1+Z1
        t13 = t12^2
        t14 = t13-YY
        Z3 = t14-ZZ 
        */

        let (x, y, z) = point;

        let x2 = x * x;
        let y2 = y * y;
        let z2 = z * z;
        let y4 = &y2 * &y2;
        let t0 = x + &y2;
        let t1 = &t0 * &t0;
        let t2 = &t1 - &x2;
        let t3 = &t2 - &y4;
        let s = felt_from_uint(2) * &t3;
        let t4 = &z2 * &z2;
        let t5 = Self::a() * &t4;
        let t6 = felt_from_uint(3) * &x2;
        let m = &t6 + &t5;
        let t7 = &m * &m;
        let t8 = felt_from_uint(2) * &s;
        let t = &t7 - &t8;
        let x3 = t.clone();
        let t9 = &s - &t;
        let t10 = felt_from_uint(8) * &y4;
        let t11 = &m * &t9;
        let y3 = &t11 - &t10;
        let t12 = y + &z;
        let t13 = &t12 * &t12;
        let t14 = &t13 - &y2;
        let z3 = &t14 - &z2;

        (x3, y3, z3)

    }
}

pub fn get_precompute() -> Vec<(FieldElement, FieldElement, FieldElement)> {
    let mut base = (G_X.clone(), G_Y.clone(), felt_from_uint(1));

    let mut precomputes = Vec::with_capacity(256);
    precomputes.push(base.clone());
    for _ in 0..256 {
        let double = Secp256k1::point_double(&base);
        precomputes.push(double.clone());
        base = double.clone();
    }

    precomputes
}

#[cfg(test)]
mod ecc_tests {
    use super::*;

    #[test]
    fn test_vectors() {
        //More to add

        let curve = Secp256k1 {};

        assert_eq!(
            curve.mul(biguint_from_str("3")),
            (
                FieldElement::new(biguint_from_str(
                    "112711660439710606056748659173929673102114977341539408544630613555209775888121"
                )),
                FieldElement::new(biguint_from_str(
                    "25583027980570883691656905877401976406448868254816295069919888960541586679410"
                ))
            )
        );
        assert_eq!(
            curve.mul(biguint_from_str("2")),
            (
                FieldElement::new(biguint_from_str(
                    "89565891926547004231252920425935692360644145829622209833684329913297188986597"
                )),
                FieldElement::new(biguint_from_str(
                    "12158399299693830322967808612713398636155367887041628176798871954788371653930"
                ))
            )
        );

        assert_eq!(curve.mul(biguint_from_str("28948022309329048855892746252171976963209391069768726095651290785379540373584")), (felt_from_str("75404758482970552478342687949548602789701733940509850780379145804275702033212"), felt_from_str("51231939447366605701190019263228486011330128519473004560491454193878655241557")));

    }

}
