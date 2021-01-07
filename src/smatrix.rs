//-------------------------------------------------------------------------
// @file smatrix.rs
//
// @date 01/06/21 14:23:23
// @author Martin Noblia
// @email mnoblia@disroot.org
//
// @brief
//
// @detail
//
//  Licence:
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or (at
// your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
//--------------------------------------------------------------------------
use std::ops::{Deref, DerefMut, Index, IndexMut};
use std::ops::{Add, Mul, Sub, AddAssign, SubAssign};

use num::{Num, Zero, One};
// use crate::traits::LinearAlgebra;

//-------------------------------------------------------------------------
//                        types
//-------------------------------------------------------------------------
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SMatrix<T, const M: usize, const N: usize>([[T; N]; M]);

impl<T, const M: usize, const N: usize> SMatrix<T, M, N> {
    pub fn new(input: [[T; N]; M]) -> Self {
        Self(input)
    }

    pub const fn shape(&self) -> (usize, usize) {
        (M, N)
    }
}

//-------------------------------------------------------------------------
//                        utils traits
//-------------------------------------------------------------------------
// indexing
impl<T, const M: usize, const N: usize> Index<(usize, usize)> for SMatrix<T, M, N> {
    type Output = T;
    fn index(&self, index: (usize, usize)) -> &T {
        &self.0[index.0][index.1]
    }
}

impl<T, const M: usize, const N: usize> IndexMut<(usize, usize)> for SMatrix<T, M, N> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut T {
        &mut self.0[index.0][index.1]
    }
}

//-------------------------------------------------------------------------
//                        rectangular matrix family M == N
//-------------------------------------------------------------------------
impl<T: Num + Copy + AddAssign, const N: usize> SMatrix<T, N, N> {
    pub fn zeros() -> Self {
        <Self as Zero>::zero()
    }

    pub fn identity() -> Self {
        <Self as One>::one()
    }
}

// Zero
impl<T: Num + Copy + AddAssign, const N: usize> Zero for SMatrix<T, N, N> {
    fn zero() -> SMatrix<T, N, N> {
        Self::new([[T::zero(); N]; N])
    }

    fn is_zero(&self) -> bool {
        *self == Self::zero()
    }
}

// One
impl<T: Num + Copy + AddAssign, const N: usize> One for SMatrix<T, N, N> {
    /// Create a identity matrix
    fn one() -> SMatrix<T, N, N> {
        let one = T::one();
        let mut result = Self::zeros();
        for i in 0..N {
            for j in 0..N {
                if i == j {
                    result[(i, j)] = one;
                }
            }
        }
        result
    }
}

// Add
impl<T: Num + Copy + AddAssign, const N: usize> Add for SMatrix<T, N, N> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self {
        let mut result = Self::zeros();
        for i in 0..N {
            for j in 0..N {
                result[(i, j)] = self[(i, j)] + rhs[(i, j)];
            }
        }
        result
    }
}

// Sub
impl<T: Num + Copy + AddAssign, const N: usize> Sub for SMatrix<T, N, N> {
    type Output = Self;
    fn sub(self, rhs: Self) -> Self {
        let mut result = Self::zeros();
        for i in 0..N {
            for j in 0..N {
                result[(i, j)] = self[(i, j)] - rhs[(i, j)];
            }
        }
        result
    }
}

// Mul
impl<T: Num + Copy + AddAssign, const N: usize> Mul for SMatrix<T, N, N> {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let mut result = SMatrix::zeros();
        for i in 0..N {
            for j in 0..N {
                for k in 0..N {
                    result[(i, j)] += self[(i, k)] * rhs[(k, j)];
                }
            }
        }
        result
    }
}

//-------------------------------------------------------------------------
//                        tests
//-------------------------------------------------------------------------
#[cfg(test)]
mod test {
    use super::SMatrix;

    #[test]
    fn add_test() {
        let x = SMatrix::new([[1, 2],
                              [1, 2]]);
        let result = x + x;
        let expected = SMatrix::new([[2, 4],
                                     [2, 4]]);
        for i in 0..result.shape().0 {
            for j in 0..result.shape().1 {
                assert_eq!(result[(i, j)], expected[(i, j)]);
            }
        }
    }

    #[test]
    fn sub_test() {
        let x = SMatrix::new([[1, 2], [1, 2]]);
        let result = x - x;
        let expected: SMatrix<u32, 2, 2> = SMatrix::zeros();
        for i in 0..result.shape().0 {
            for j in 0..result.shape().1 {
                assert_eq!(result[(i, j)], expected[(i, j)]);
            }
        }
    }

    #[test]
    fn zeros_test() {
        let x:SMatrix<u32,100,100> = SMatrix::zeros();
        for i in 0..x.shape().0 {
            for j in 0..x.shape().1 {
                assert_eq!(x[(i, j)], 0);
            }
        }
    }

    #[test]
    fn mul_test() {
        let x = SMatrix::new([[1, 2], [1, 2]]);
        let result = x * x;
        let expected = SMatrix::new([[3, 6], [3, 6]]);
        for i in 0..x.shape().0 {
            for j in 0..x.shape().1 {
                assert_eq!(result[(i, j)], expected[(i, j)]);
            }
        }
    }
}

