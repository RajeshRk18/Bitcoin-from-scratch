use num_bigint::BigInt;
use crate::ecc::PointError;

pub trait EllipticCurve {
    fn is_on_curve(&self, point: (BigInt, BigInt)) -> bool;

    fn check_point_not_at_infinity(
        &self,
        point1: (BigInt, BigInt),
        point2: (BigInt, BigInt),
    ) -> bool;

    fn point_add(
        &self,
        point1: (BigInt, BigInt),
        point2: (BigInt, BigInt),
    ) -> Result<(BigInt, BigInt), PointError>;

    fn point_double(&self, point: (BigInt, BigInt)) -> (BigInt, BigInt);

    fn scalar_mul(&self, scalar: BigInt, point: (BigInt, BigInt)) -> Result<(BigInt, BigInt), PointError>;
}

   /*impl TestCurve {
        fn ext_euclid(a: BigInt, b: BigInt) -> (BigInt, BigInt, BigInt) {
            if b == 0i8.to_bigint().unwrap() {
                return (a.clone(), 1i8.to_bigint().unwrap(), 0i8.to_bigint().unwrap());
            } else {
                let (gcd, x1, y1) = Self::ext_euclid(b.clone(), a.clone() % b.clone());
                let x = y1.clone();
                let y = x1.clone() - (a / b) * y1.clone();
                return (gcd, x, y);
            }
        }
    
        fn mod_inv(&self, a: BigInt, b: BigInt) -> Option<BigInt> {
            let (gcd, x, _) = Self::ext_euclid(a, b.clone());
    
            if gcd != 1.to_bigint().unwrap() {
                return None;
            } else {
                return Some((x % b.clone() + b.clone()) % b.clone());
            }
        }
    }

    impl EllipticCurve for TestCurve {
        fn is_on_curve(&self, point: (BigInt, BigInt)) -> bool {
            let lhs = point.1.pow(2) % BigInt::from_i128(self.order).unwrap();
            let rhs = (point.0.pow(3) + self.b) % BigInt::from_i128(self.order).unwrap();
            if lhs != rhs {
                return false;
            }
            true
        }

        fn check_point_not_at_infinity(
            &self,
            point1: (BigInt, BigInt),
            point2: (BigInt, BigInt),
        ) -> bool {
            if point1.0 == point2.0 {
                return false;
            }
            true
        }

        fn point_add(
            &self,
            point1: (BigInt, BigInt),
            point2: (BigInt, BigInt),
        ) -> Result<(BigInt, BigInt), PointError> {
            if !self.is_on_curve(point1.clone()) {
                return Err(PointError(point2.0.clone(), point2.1.clone()));
            }
            if !self.is_on_curve(point2.clone()) {
                return Err(PointError(point2.0.clone(), point2.1.clone()));
            }
            let (x1, y1) = point1;
            let (x2, y2) = point2;

            if x1 != x2 {
                let slope = (y2.clone() - y1.clone()) / (x2.clone() - x1.clone());
                let x3 =
                    (slope.pow(2) - x1.clone() - x2.clone()) % BigInt::from(self.order as isize);
                let y3 = (slope * (x1.clone() - x3.clone()) - y1.clone())
                    % BigInt::from(self.order as isize);
                return Ok((x3, y3));
            } else if ((x1 == x2) && (y1 == -y2.clone())) || ((x1 == x2) && (y1 != y2)) {
                return Ok((BigInt::from(1 as usize), BigInt::from(0 as usize)));
            } else if ((x1 == x2) && (y1 == y2)) || ((y1 == y2) && (x1 != x2)) {
                let slope: BigInt = (BigInt::from_u32(3).unwrap() * x1.clone().pow(2)
                    / BigInt::from_u32(2).unwrap()
                    * y1.clone())
                    % BigInt::from(self.order as isize);
                let x3: BigInt = (slope.clone().pow(2) - BigInt::from_u32(2).unwrap() * x1.clone())
                    % BigInt::from(self.order as isize);
                let y3: BigInt = (slope * (x1.clone() - x3.clone()) - y1.clone())
                    % BigInt::from(self.order as isize);
                return Ok((x3, y3));
            } else {
                return Ok((BigInt::from(1 as usize), BigInt::from(0 as usize)));
            }
        }

        fn point_double(&self, point: (BigInt, BigInt)) -> (BigInt, BigInt) {
            if point.1 == Self::zero() {
                return Self::point_at_infinity();
            } else {
                let slope = (BigInt::from_u8(3).unwrap() * point.0.clone() + self.degree0_coeff.clone()) / BigInt::from_u8(2).unwrap() * point.1.clone();
                let x3 = slope.pow(2) - BigInt::from_u8(2).unwrap() * point.0.clone();
                let y3 = point.1 + slope * (x3.clone() - point.0);
                return (x3, y3);
            }
        }

        fn scalar_mul(&self, mut scalar: BigInt, point: (BigInt, BigInt)) -> Result<(BigInt, BigInt), PointError> {
            assert!(scalar.clone() != Self::zero());
            let mut q = (BigInt::from_u8(0).unwrap(), BigInt::from_u8(0).unwrap());
    
            let mut r = point.clone();
    
            while scalar.clone() > BigInt::from_u8(0).unwrap() {
                if scalar.clone() % 2 == Self::one() {
                    q = self.point_add(q.clone(), r.clone())?;
                } else {
                    scalar = scalar.clone() >> 1;
                    r = self.point_double(r.clone());
                }
            }
    
            Ok(q)
        }

    }

    #[test]
    fn in_curve() {
        let test_curve = TestCurve { b: 7, order: 223 };
        assert!(test_curve.is_on_curve((
            BigInt::from_u32(192).unwrap(),
            BigInt::from_u32(105).unwrap()
        )));
        assert!(!test_curve.is_on_curve((
            BigInt::from_u32(200).unwrap(),
            BigInt::from_u32(119).unwrap()
        )));
    }

    #[test]
    fn add() {
        let test_curve = TestCurve { b: 7, order: 223 };
        debug_assert_eq!(
            test_curve.point_add(
                (
                    BigInt::from_u32(192).unwrap(),
                    BigInt::from_u32(105).unwrap()
                ),
                (BigInt::from_u32(17).unwrap(), BigInt::from_u32(56).unwrap())
            ),
            Ok((
                BigInt::from_i32(-209).unwrap(),
                BigInt::from_i32(-105).unwrap()
            ))
        );
    }

   /*#[test]
    fn mul() {
        let test_curve = TestCurve { b: 7, order: 223 };
        debug_assert_eq!(
            test_curve.scalar_mul(
                3,
                (
                    BigInt::from_u32(192).unwrap(),
                    BigInt::from_u32(105).unwrap()
                )
            ),
            (
                BigInt::from_i32(18).unwrap(),
                BigInt::from_i32(189).unwrap()
            )
        );
    }*/*/
