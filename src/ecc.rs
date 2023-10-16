use crate::fieldelement::*;
use crate::traits::{IsEllipticCurve, Secp256k1Curve, Signable};
use crate::utils::*;
use lazy_static::lazy_static;
use num_bigint::BigUint;
use num_traits::FromPrimitive;
use rand_chacha::ChaCha20Rng;

lazy_static! {
    pub static ref BASE_ORDER: FieldElement = FieldElement::new(
        BigUint::parse_bytes(
            "115792089237316195423570985008687907852837564279074904382605163141518161494337"
                .as_bytes(),
            10
        )
        .unwrap()
    );
    pub static ref CURVE_ORDER: FieldElement = FieldElement::new(
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

    pub fn generator_proj() -> <Self as IsEllipticCurve>::Jacobian {
        (G_X.clone(), G_Y.clone(), felt_from_uint(1))
    }

    pub fn generator_affine() -> <Self as IsEllipticCurve>::Affine {
        Self::generator()
    }

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
        BASE_ORDER.clone()
    }

    fn scalar_field_order() -> FieldElement {
        CURVE_ORDER.clone()
    }

    fn generator() -> Self::Affine {
        (G_X.clone(), G_Y.clone())
    }

    fn identity() -> Self::Jacobian {
        (FieldElement::one(), FieldElement::one(), FieldElement::zero())
    }
}

impl IsEllipticCurve for Secp256k1 {
    type Affine = (FieldElement, FieldElement);
    type Jacobian = (FieldElement, FieldElement, FieldElement);
    type Scalar = BigUint;

    #[inline(always)]
    fn is_valid(&self, point: &Self::Affine) -> bool {
        let lhs = &point.1 * &point.1;
        let rhs = &point.0 * &point.0 * &point.0 + felt_from_uint(7);
        
        lhs == rhs
    }

    // secp256k1 params have been chosen to be non-discriminant so that it is smooth and non-singular. 
    fn is_indiscriminant() -> bool {
        let op = -felt_from_uint(16)
            * (felt_from_uint(4) * Self::a() * Self::a() * Self::a()
                + felt_from_uint(27) * Self::b() * Self::b());

        op != FieldElement::zero()
    }

    fn is_negative(point1: &Self::Jacobian, point2: &Self::Jacobian) -> bool {
        (&point1.0 == &point2.0) && (&point1.2 == &point2.2) && (&point1.1 == &(-&point2.1)) 
    }

    // checks if the point is additive identity.
    fn is_identity(p: &Self::Jacobian) -> bool {
        p.2 == FieldElement::zero()
    }

    // Performs Group Isomorphism: Maps point in Affine Space to Projective.
    // (x, y) -> (x, y, 1)
    #[inline(always)]
    fn to_jacobian(&self, point: Self::Affine) -> Self::Jacobian {
        (point.0, point.1, FieldElement::one())
    }

    // Maps point in Projective Space to Affine.
    // (x, y, z) -> (x, y)
    #[inline(always)]
    fn to_affine(&self, point: Self::Jacobian) -> Self::Affine {
        let (x, y, z) = point;

        let z_inv = &FieldElement::one() / &z;
        let z_inv_squared = &z_inv * &z_inv;
        let z_inv_cubed = &z_inv_squared * &z_inv;
        let x = x * z_inv_squared;
        let y = y * z_inv_cubed;

        if z == FieldElement::zero() {
            return (FieldElement::zero(), FieldElement::zero());
        } else {
            return (x, y);
        }
    }

    #[inline(always)]
    fn scalar_gen_mul(&self, mut scalar: BigUint) -> Self::Jacobian {
        assert!(
            scalar != biguint_from_uint(0 as u32)
                && FieldElement::from(scalar.clone()) < BASE_ORDER.clone()
        );
        let mut addend = Self::identity();
        //let precomputes = Self::get_precomputed_table();

        let mut result: (FieldElement, FieldElement, FieldElement) = Self::generator_proj();
        let mut dummy_point: (FieldElement, FieldElement, FieldElement) = Self::identity();

        for _ in 0 .. scalar.bits() {
            if scalar.clone() & BigUint::from_u8(1).unwrap() == BigUint::from_u8(1).unwrap() {
                addend = Self::point_add(&addend, &result);
                dbg!(&result);
            } else {
                // Dummy op to make it constant time.
                dummy_point = Self::point_add(&dummy_point, &dummy_point);
            }
            result = Self::point_double(&result);
            dbg!(self.to_affine(result.clone()));
            scalar >>= 1;
        }
        addend
    }

    fn point_add(point1: &Self::Jacobian, point2: &Self::Jacobian) -> Self::Jacobian {
        let mut t = Vec::with_capacity(5);

        for _ in 0..5 {
            t.push(Self::identity());
        }
        
        t[1] = point2.clone();
        t[4] = point1.clone();

        let p1_is_identity = Self::is_identity(&point1);
        let p2_is_identity = Self::is_identity(&point2);
        
        let (x1, y1, z1) = point1;
        let (x2, y2, z2) = point2;

        let mut t2 = z1 * z1;
        let mut t3 = z1 * &t2;
        let mut t1 = x2 * &t2;
        let mut t4 = y2 * &t3;
        t3 = z2 * z2;
        let mut t5 = z2 * &t3;
        let mut t7 = x1 * &t3;
        let mut t8 = y1 * &t5;
        let mut t1 = t1 - &t7;
        t4 -= t8.clone();

        let mut i = 3;

        if t1 == FieldElement::zero() {
            i = 0;
            
            if t4 == FieldElement::zero() {
                i = 2;
            }
        }

        if p1_is_identity {
            i = 1;
        }

        if p2_is_identity {
            i = 4;
        }

        let mut mask = 0;

        if i == 3 {
            mask = 1;
        }

        let mut t3 = x1 + &t1;
        let mut t6 = x1 - &t2;

        if mask == 0 {
            t2 = y1.clone();
        } else {
            t2 = t1.clone();
        }

        let t5 = &t2 * &t2;
        
        if mask == 0 {
            t7 = x1.clone();
        }

        t1 = &t5 * &t7;

        let table_z = z1 * &t2;
        t[2].2 = table_z;

        let table_z = z2 * z2;
        t[3].2 = table_z;

        if mask != 0 {
            t3 = t2.clone();
            t6 = t5.clone();
        }

        t2 = &t3 * &t6;

        t3 = FieldElement::floor_div(&t2, &felt_from_uint(2));
        t3 += t2.clone();

        if mask != 0 {
            t3 = t4.clone();
        }

        t4 = &t3 * &t3;
        t4 -= t1.clone();

        t[2].0 = &t4 - &t1;
        t[3].0 = &t[2].0 - &t2;

        if mask == 0 {
            t4 = t[2].0.clone();
        } else {
            t4 = t[3].0.clone();
        }

        t1 -= t4.clone();
        t4 = &t3 * &t1;

        if mask == 0 {
            t1 = t5.clone();
            t2 = t5.clone();
        } else {
            t1 = t8.clone();
        }

        t3 = &t1 * &t2;
        t[2].1 = &t4 - &t3;
        t[3].1 = t[2].1.clone();

        t[i].clone()


            /*if Self::is_identity(&point2) {
                point1.clone()
            } else if Self::is_identity(&point1) {
                point2.clone()
            } else {
                let (px, py, pz) = point1;
                let (qx, qy, qz) = point2;

                let z1z1 = pz * pz;
                let u2 = qx * pz;
                let s2 = qy * pz;
                let s2 = s2 * &z1z1;

                if &u2 == pz && &s2 == py {
                    return Self::point_double(point1);
                }

                let h = &u2 - px;
                let hh = &h * &h;
                let i = &hh + &hh;
                let i = &i + &i;
                let j = &h * &i;
                let r = &s2 - &py;
                let r = &r * &r;
                let v = px * &i;
                let px = &r * &r;
                let px = px - &j;
                let px = px - &v;
                let j = j * py;
                let j = &j * &j;
                let py = &v - &px;
                let py = py * &r;
                let py = py - &j;
                let pz = pz + &h;
                let pz = &pz * &pz;
                let pz = pz - &z1z1;
                let pz = pz - hh;

                return (px, py, pz);*/



                /*let u1 = qy * pz;
                let u2 = py * qz;
                let v1 = qx * pz;
                let v2 = px * qz;
                if v1 == v2 {
                    if u1 != u2 || *py == FieldElement::zero() {
                        Self::identity()
                    } else {
                        let px_square = px * px;
                        let three_px_square = &px_square + &px_square + &px_square;
                        let w = Self::a() * pz * pz + three_px_square;
                        let w_square = &w * &w;
    
                        let s = py * pz;
                        let s_square = &s * &s;
                        let s_cube = &s * &s_square;
                        let two_s_cube = &s_cube + &s_cube;
                        let four_s_cube = &two_s_cube + &two_s_cube;
                        let eight_s_cube = &four_s_cube + &four_s_cube;
    
                        let b = px * py * &s;
                        let two_b = &b + &b;
                        let four_b = &two_b + &two_b;
                        let eight_b = &four_b + &four_b;
    
                        let h = &w_square - &eight_b;
                        let hs = &h * &s;
    
                        let pys_square = py * py * s_square;
                        let two_pys_square = &pys_square + &pys_square;
                        let four_pys_square = &two_pys_square + &two_pys_square;
                        let eight_pys_square = &four_pys_square + &four_pys_square;
    
                        let xp = &hs + &hs;
                        let yp = w * (four_b - &h) - eight_pys_square;
                        let zp = eight_s_cube;
                        (xp, yp, zp)
                    }
                } else {
                    let u = u1 - &u2;
                    let v = v1 - &v2;
                    let w = pz * qz;
    
                    let u_square = &u * &u;
                    let v_square = &v * &v;
                    let v_cube = &v * &v_square;
                    let v_square_v2 = &v_square * &v2;
    
                    let a = &u_square * &w - &v_cube - (&v_square_v2 + &v_square_v2);
    
                    let xp = &v * &a;
                    let yp = u * (&v_square_v2 - &a) - &v_cube * &u2;
                    let zp = &v_cube * &w;
                    (xp, yp, zp)
                }*/
            }

    /*#[inline(always)]
    fn point_add(point1: &Self::Jacobian, point2: &Self::Jacobian) -> Self::Jacobian {

    /*
     Z1Z1:=Z1^2;
     U2:=X2*Z1Z1;
     S2:=Y2*Z1*Z1Z1;
     H:=U2-X1;
     HH:=H^2;
     I:=4*HH;
     J:=H*I;
     r:=2*(S2-Y1);
     V:=X1*I;
     X3:=r^2-J-2*V;
     Y3:=r*(V-X3)-2*Y1*J;
     Z3:=(Z1+H)^2-Z1Z1-HH;
    */
        //dbg!(&point1, point2);
        let is_negative = Self::is_negative(&point1, &point2);
        dbg!(is_negative);
        #[allow(non_snake_case)]
        let P = point1.clone();
        
        let dummy_val = Self::identity();

        let (x1, y1, z1) = point1;
        let (x2, y2, z2) = point2;

        let is_equal = x1 == x2 && y1 == y2 && z1 == z2;

        let z1z1 = z1 * z1;
        dbg!(&z1z1);
        let u = x2 * &z1z1;
        dbg!(&u);
        let s = y2 * z1 * &z1z1;
        dbg!(&s);
        let h = &u - x1;
        dbg!(&h);
        let hh = &h * &h;
        dbg!(&hh);
        let i = felt_from_uint(4) * &hh;
        dbg!(&i);
        let j = &h * &i;
        dbg!(&j);
        let r = felt_from_uint(2) * (&s - &y1);
        dbg!(&r);
        let v = x1 * &i;
        dbg!(&v);
        let x3 = &r * &r - &j - felt_from_uint(2) * &v;
        dbg!(&x3);
        let y3 = &r * &(&v - &x3) - felt_from_uint(2) * y1 * j;
        dbg!(&y3);
        let z3 = (z1 + &h) * (z1 + &h) - z1z1 + hh;
        dbg!(&z3);
        
        if is_equal {
            return Self::point_double(&P);
        } else if is_negative {
            #[allow(unused_variables)]
            let dummy_op = Self::point_double(&dummy_val);
            return Self::identity();
        } else {
            #[allow(unused_variables)]
            let dummy_op = Self::point_double(&dummy_val);    
            return (x3, y3, z3);        
        }
    }*/

    #[inline(always)]
    fn point_double(point: &Self::Jacobian) -> Self::Jacobian {
        // https://www.hyperelliptic.org/EFD/oldefd/jacobian.html#DBL

        /*
             XX:=X1^2;
     YY:=Y1^2;
     ZZ:=Z1^2;
     S:=4*X1*YY;
     M:=3*XX+a*ZZ^2;
     T:=M^2-2*S;
     X3:=T;
     Y3:=M*(S-T)-8*YY^2;
     Z3:=2*Y1*Z1; 
     */

        /*let (x, y, z) = point;

        let zz = z * z;
        let 
        let yy = y * y;
        let yyyy = &yy * &yy;
        let xx = x * x;

        let s = x + &yy;
        let s = (&s * &s - &xx - &yyyy) * felt_from_uint(2);
        let m = &xx * &xx + &xx;
        let z = (z + y) * (z + y) - &yy - &zz;
        let t = &m * &m;
        let x = t.clone();
        let t = &s * &s;
        let x = &x - &t;
        let y = (s - &x) * &m;
        let yyyy = &yyyy * &yyyy;
        let yyyy = &yyyy * &yyyy;
        let yyyy = &yyyy * &yyyy;
        let y = y - yyyy;

        return (x, y, z);
*/
        let (x, y, z) = point;

        let t1 = z * z;
        let t2 = x + &t1;
        let t1 = x - &t1;
        let t1  = t1 * t2;
        let t2 =  FieldElement::floor_div(&t1, &felt_from_uint(2));
        let t1 = &t1 + &t2;
        let t2 = y * y;
        let t3 = x * &t2;
        let t4 = &t1 * &t1;
        let t4 = t4 - &t3;
        let x1 = &t4 - &t3;
        let z1 = y * z;
        let t2 = &t2 * &t2; 
        let t4 = t3 - &x1;
        let t1 = t1 * t4;
        let y1 = t1 - t2;

        (x1, y1, z1)

        
        /*dbg!(&point);        
        let (x, y, z) = point;
        dbg!(&x);
        let xx = x * x;
        dbg!(&xx);
        let yy = y * y;
        dbg!(&yy);
        let s = felt_from_uint(4) * x * &yy;
        dbg!(&s);
        let m = felt_from_uint(3) * xx;
        dbg!(&m);
        let t = &m * &m - felt_from_uint(2) * &s;
        dbg!(&t);
        let x2 = t.clone();
        dbg!(&x2);
        let y2 = m * (s - t) - felt_from_uint(8) * &yy * &yy;
        dbg!(&y2);
        let z2 = felt_from_uint(2) * y * z;
        dbg!(&z2);

        (x2, y2, z2)*/
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
mod ecc_test {
    use super::*;

    #[test]
    fn test_vectors() {
        //More to add

        let curve = Secp256k1 {};

        assert_eq!(curve.to_affine(Secp256k1::point_double(&Secp256k1::generator_proj())), Secp256k1::generator());

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

        assert_eq!(
            curve.mul(biguint_from_str(
                "28948022309329048855892746252171976963209391069768726095651290785379540373584"
            )),
            (
                felt_from_str(
                    "75404758482970552478342687949548602789701733940509850780379145804275702033212"
                ),
                felt_from_str(
                    "51231939447366605701190019263228486011330128519473004560491454193878655241557"
                )
            )
        );
    }
}
