use crate::fieldelement::*;
use num_bigint::{BigUint, ToBigUint};

pub fn felt_from_uint<T: ToBigUint>(value: T) -> FieldElement {
    FieldElement::new(value.to_biguint().unwrap())
}

pub fn felt_from_str(value: &str) -> FieldElement {
    FieldElement::new(BigUint::parse_bytes(value.as_bytes(), 10).unwrap())
}

// As we work under prime field, there always exists multiplicative inverse. So, no need to wrap output with Option.
pub fn ext_euclid(a: &FieldElement, b: &FieldElement) -> FieldElement {
    let mut a = a.clone();
    let mut b = b.clone();
    if gcd(&a, &b) != 1 {
        return felt_from_uint(0);
    } else {
        let mut s = felt_from_uint(0);
        let mut s_old = felt_from_uint(1);
        let mut t = felt_from_uint(1);
        let mut t_old = felt_from_uint(0);

        while b != 0 {
            // Floor div
            let q = FieldElement::floor_div(&a, &b);

            let temp = b.clone();
            b = a - &q * &b;
            a = temp;

            let temp = s.clone();
            s = s_old - &q * &s;
            s_old = temp;

            let temp = t.clone();
            t = t_old - &q * &t;
            t_old = temp;
        }

        s_old
    }
}

pub fn gcd(a: &FieldElement, b: &FieldElement) -> FieldElement {
    let mut a = a.clone();
    let mut b = b.clone();
    while b != 0 {
        let t = b.clone();
        b = a % b;
        a = t;
    }
    return a;
}

#[allow(non_snake_case)]
pub fn biguint_from_uint<T: ToBigUint>(num: T) -> BigUint
where
    BigUint: From<T>,
{
    num.into()
}

#[allow(non_snake_case)]
pub fn biguint_from_str(x: &str) -> BigUint {
    BigUint::parse_bytes(x.as_bytes(), 10).expect("cannot create Biguint from the gn input")
}

#[cfg(test)]
mod util_test {
    use super::*;
    use crate::ecc::CURVE_ORDER;

    #[test]
    fn inv_test() {
        let a = felt_from_str(
            "18233444644265725414720095600783733811551140822142748479967321742263849032310",
        );
        debug_assert_eq!(
            ext_euclid(&a, &CURVE_ORDER.clone()),
            felt_from_str(
                "7137376184023522026654217343440040540351594155614695530712440441403759244341"
            )
        );
        debug_assert_eq!(
            ext_euclid(&felt_from_uint(6547445), &felt_from_str("7855654643")),
            felt_from_str("3482694415")
        );
    }
}
