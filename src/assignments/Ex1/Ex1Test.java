package assignments.Ex1;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

/**
 *  * Introduction to Computer Science 2026, Ariel University,
 *  * Ex1: arrays, static functions and JUnit
 *
 * This JUnit class represents a JUnit (unit testing) for Ex1-
 * It contains few testing functions for the polynomial functions as define in Ex1.
 * Note: you should add additional JUnit testing functions to this class.
 *
 * @author boaz.ben-moshe
 */

class Ex1Test {
	static final double[] P1 ={2,0,3, -1,0}, P2 = {0.1,0,1, 0.1,3};
	static double[] po1 = {2,2}, po2 = {-3, 0.61, 0.2};;
	static double[] po3 = {2,1,-0.7, -0.02,0.02};
	static double[] po4 = {-3, 0.61, 0.2};
	
 	@Test
	/**
	 * Tests that f(x) == poly(x).
	 */
	void testF() {
		double fx0 = Ex1.f(po1, 0);
		double fx1 = Ex1.f(po1, 1);
		double fx2 = Ex1.f(po1, 2);
		assertEquals(fx0, 2, Ex1.EPS);
		assertEquals(fx1, 4, Ex1.EPS);
		assertEquals(fx2, 6, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1(x) + p2(x) == (p1+p2)(x)
	 */
	void testF2() {
		double x = Math.PI;
		double[] po12 = Ex1.add(po1, po2);
		double f1x = Ex1.f(po1, x);
		double f2x = Ex1.f(po2, x);
		double f12x = Ex1.f(po12, x);
		assertEquals(f1x + f2x, f12x, Ex1.EPS);
	}
	@Test
	/**
	 * Tests that p1+p2+ (-1*p2) == p1
	 */
	void testAdd() {
		double[] p12 = Ex1.add(po1, po2);
		double[] minus1 = {-1};
		double[] pp2 = Ex1.mul(po2, minus1);
		double[] p1 = Ex1.add(p12, pp2);
		assertTrue(Ex1.equals(p1, po1));
	}
	@Test
	/**
	 * Tests that p1+p2 == p2+p1
	 */
	void testAdd2() {
		double[] p12 = Ex1.add(po1, po2);
		double[] p21 = Ex1.add(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
	@Test
	/**
	 * Tests that p1+0 == p1
	 */
	void testAdd3() {
		double[] p1 = Ex1.add(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, po1));
	}

    //-------------------------tests add--------------------------------
    @Test
    void testCompactRemovesTrailingZeros() {
        double[] p = {1, 2, 3, 0, 0};           /// [1] Polynomial with unnecessary zeros at the end
        double[] expected = {1, 2, 3};          /// [2] Expected result after compact: only non-zero coefficients remain
        assertTrue(Ex1.equals(Ex1.compact(p), expected)); /// [3] Verify compact removes trailing zeros correctly
    }


    @Test
    void testCompactAllZeros() {
        double[] p = {0, 0, 0};                 /// [1] Polynomial with only zeros
        double[] expected = {0};                /// [2] Expected result: single zero (representing zero polynomial)
        assertTrue(Ex1.equals(Ex1.compact(p), expected)); /// [3] Verify compact handles all-zero input correctly
    }


    @Test
    void testAddSameLengthAfterCompact() {
        double[] p1 = {1, 2, 3, 0, 0};          /// [1] Polynomial with trailing zeros
        double[] p2 = {4, 5, 6, 0};             /// [2] Another polynomial with trailing zeros
        double[] expected = {5, 7, 9};          /// [3] Expected sum after compact: (1+4), (2+5), (3+6)
        assertTrue(Ex1.equals(Ex1.add(p1, p2), expected)); /// [4] Verify add works correctly after compacting both inputs
    }


    @Test
    void testAddDifferentLengthAfterCompact() {
        double[] p1 = {1, 2, 3, 4, 0, 0};       /// [1] Longer polynomial with trailing zeros
        double[] p2 = {5, 6, 0};                /// [2] Shorter polynomial with trailing zeros
        double[] expected = {1, 2, 8, 10};      /// [3] Expected sum: first two terms stay, then (3+5)=8, (4+6)=10
        assertTrue(Ex1.equals(Ex1.add(p1, p2), expected)); /// [4] Verify add handles different lengths after compact
    }


    @Test
    void testAddWithNullInputs() {
        double[] p1 = null;                     /// [1] First polynomial is null
        double[] p2 = {1, 2, 0};                /// [2] Second polynomial with trailing zero
        double[] expected = {1, 2};             /// [3] Expected result after compact: {1, 2}
        assertTrue(Ex1.equals(Ex1.add(p1, p2), expected)); /// [4] Verify add returns p2 when p1 is null
    }

///
    @Test
    void testAddBothNull() {
        double[] p1 = null;                     /// [1] First polynomial is null
        double[] p2 = null;                     /// [2] Second polynomial is null
        double[] expected = {0};                /// [3] Expected result: ZERO polynomial
        assertTrue(Ex1.equals(Ex1.add(p1, p2), expected)); /// [4] Verify add returns ZERO when both inputs are null
    }


    @Test
    void testAddEmptyArrays() {
        double[] p1 = {};                       /// [1] First polynomial is empty
        double[] p2 = {1, 2, 0};                /// [2] Second polynomial with trailing zero
        double[] expected = {1, 2};             /// [3] Expected result after compact: {1, 2}
        assertTrue(Ex1.equals(Ex1.add(p1, p2), expected)); /// [4] Verify empty array behaves like zero
    }


    @Test
    void testAddNegativeValuesAfterCompact() {
        double[] p1 = {-1, -2, -3, 0};          /// [1] Polynomial with negative coefficients and trailing zero
        double[] p2 = {4, 5, 6, 0};             /// [2] Polynomial with positive coefficients and trailing zero
        double[] expected = {3, 3, 3};          /// [3] Expected sum: (-1+4), (-2+5), (-3+6)
        assertTrue(Ex1.equals(Ex1.add(p1, p2), expected)); /// [4] Verify add works with negatives after compact
    }
//-------------------------end of tests add--------------------------------

    @Test
	/**
	 * Tests that p1*0 == 0
	 */
	void testMul1() {
		double[] p1 = Ex1.mul(po1, Ex1.ZERO);
		assertTrue(Ex1.equals(p1, Ex1.ZERO));
	}
	@Test
	/**
	 * Tests that p1*p2 == p2*p1
	 */
	void testMul2() {
		double[] p12 = Ex1.mul(po1, po2);
		double[] p21 = Ex1.mul(po2, po1);
		assertTrue(Ex1.equals(p12, p21));
	}
    /// -------------------------tests mul--------------------------------
    ///
    @Test
    void testMulBasic() {
        double[] p1 = {1, 2};                     // [1] Polynomial: 2x + 1
        double[] p2 = {3, 4};                     // [2] Polynomial: 4x + 3
        double[] expected = {3, 10, 8};           // [3] Result: (1*3), (1*4 + 2*3), (2*4) → 3 + 10x + 8x²
        assertTrue(Ex1.equals(Ex1.mul(p1, p2), expected)); // [4] Verify multiplication works for basic case
    }


    @Test
    void testMulWithNull() {
        double[] p1 = null;                       // [1] First polynomial is null
        double[] p2 = {1, 2};                     // [2] Second polynomial: 2x + 1
        double[] expected = null;                 // [3] Multiplication with null should return null
        assertNull(Ex1.mul(p1, p2));              // [4] Verify result is null
    }


    @Test
    void testMulBothNull() {
        double[] p1 = null;                       // [1] First polynomial is null
        double[] p2 = null;                       // [2] Second polynomial is null
        double[] expected = null;                 // [3] Multiplication with both null should return null
        assertNull(Ex1.mul(p1, p2));              // [4] Verify result is null
    }


    @Test
    void testMulWithZeroPolynomial() {
        double[] p1 = {0};                        // [1] First polynomial is ZERO
        double[] p2 = {1, 2, 3};                  // [2] Second polynomial: 3x² + 2x + 1
        double[] expected = {0};                  // [3] ZERO times anything is ZERO
        assertTrue(Ex1.equals(Ex1.mul(p1, p2), expected)); // [4] Verify multiplication with ZERO returns ZERO
    }


    @Test
    void testMulWithIdentityPolynomial() {
        double[] p1 = {1};                        // [1] First polynomial is identity: 1
        double[] p2 = {4, 5, 6};                  // [2] Second polynomial: 6x² + 5x + 4
        double[] expected = {4, 5, 6};            // [3] 1 * p2 = p2
        assertTrue(Ex1.equals(Ex1.mul(p1, p2), expected)); // [4] Verify multiplication with identity returns other polynomial
    }


    @Test
    void testMulWithTrailingZeros() {
        double[] p1 = {1, 2, 0};                  // [1] Polynomial with trailing zero
        double[] p2 = {3, 0};                     // [2] Polynomial with trailing zero
        double[] expected = {3, 6, 0};            // [3] Result after compact: (1*3), (1*0 + 2*3), (2*0)
        assertTrue(Ex1.equals(Ex1.mul(p1, p2), expected)); // [4] Verify compact works before multiplication
    }


    @Test
    void testMulNegativeCoefficients() {
        double[] p1 = {-1, 2};                    // [1] Polynomial: 2x - 1
        double[] p2 = {3, -4};                    // [2] Polynomial: -4x + 3
        double[] expected = {-3, 10, -8};         // [3] Result: (-1*3), (-1*-4 + 2*3), (2*-4)
        assertTrue(Ex1.equals(Ex1.mul(p1, p2), expected)); // [4] Verify multiplication works with negatives
    }


    @Test
    void testMulEmptyArrays() {
        double[] p1 = {};                         // [1] First polynomial is empty
        double[] p2 = {1, 2};                     // [2] Second polynomial: 2x + 1
        double[] expected = null;                 // [3] Multiplication with empty array should return null
        assertNull(Ex1.mul(p1, p2));              // [4] Verify result is null
    }


    @Test
    void testMulLargePolynomials() {
        double[] p1 = {1, 2, 3};                  // [1] Polynomial: 3x² + 2x + 1
        double[] p2 = {4, 5};                     // [2] Polynomial: 5x + 4
        double[] expected = {4, 13, 22, 15};      // [3] Result: multiply term by term
        assertTrue(Ex1.equals(Ex1.mul(p1, p2), expected)); // [4] Verify multiplication works for larger polynomials
    }

    /// -------------------------end of tests mul--------------------------------
	@Test
	/**
	 * Tests that p1(x) * p2(x) = (p1*p2)(x),
	 */
	void testMulDoubleArrayDoubleArray() {
		double[] xx = {0,1,2,3,4.1,-15.2222};
		double[] p12 = Ex1.mul(po1, po2);
		for(int i = 0;i<xx.length;i=i+1) {
			double x = xx[i];
			double f1x = Ex1.f(po1, x);
			double f2x = Ex1.f(po2, x);
			double f12x = Ex1.f(p12, x);
			assertEquals(f12x, f1x*f2x, Ex1.EPS);
		}
	}

	@Test
	/**
	 * Tests a simple derivative examples - till ZERO.
	 */
	void testDerivativeArrayDoubleArray() {
		double[] p = {1,2,3}; // 3X^2+2x+1
		double[] pt = {2,6}; // 6x+2
		double[] dp1 = Ex1.derivative(p); // 2x + 6
		double[] dp2 = Ex1.derivative(dp1); // 2
		double[] dp3 = Ex1.derivative(dp2); // 0
		double[] dp4 = Ex1.derivative(dp3); // 0
		assertTrue(Ex1.equals(dp1, pt));
		assertTrue(Ex1.equals(Ex1.ZERO, dp3));
		assertTrue(Ex1.equals(dp4, dp3));
	}
    //-------------------------tests derivative--------------------------------

    @Test
    void testDerivativeNullInput() {
        double[] p = null;                           /// [1] Polynomial is null
        double[] expected = {0};                    /// [2] Expected result: ZERO polynomial
        assertTrue(Ex1.equals(Ex1.derivative(p), expected)); /// [3] Verify derivative of null returns ZERO
    }

    @Test
    void testDerivativeConstantPolynomial() {
        double[] p = {5};                           /// [1] Polynomial is constant: 5
        double[] expected = {0};                    /// [2] Derivative of constant is zero
        assertTrue(Ex1.equals(Ex1.derivative(p), expected)); /// [3] Verify derivative of constant returns ZERO
    }

    @Test
    void testDerivativeEmptyArray() {
        double[] p = {};                            // [1] Polynomial is empty
        double[] expected = {0};                    /// [2] Expected result: ZERO polynomial
        assertTrue(Ex1.equals(Ex1.derivative(p), expected)); /// [3] Verify derivative of empty array returns ZERO
    }

    @Test
    void testDerivativeWithZeros() {
        double[] p = {1, 0, 3};                     /// [1] Polynomial: 3x² + 0x + 1
        double[] expected = {0, 6};                 /// [2] Derivative: 6x + 0
        assertTrue(Ex1.equals(Ex1.derivative(p), expected)); /// [3] Verify derivative handles zeros correctly
    }

    @Test
    void testDerivativeNegativeCoefficients() {
        double[] p = {-1, -2, -3};                  /// [1] Polynomial: -3x² - 2x - 1
        double[] expected = {-2, -6};               /// [2] Derivative: -6x - 2
        assertTrue(Ex1.equals(Ex1.derivative(p), expected)); /// [3] Verify derivative works with negative coefficients
    }


    @Test
    void testDerivativeLargePolynomial() {
        double[] p = {1, 2, 3, 4, 5};               /// [1] Polynomial: 5x⁴ + 4x³ + 3x² + 2x + 1
        double[] expected = {2, 6, 12, 20};         /// [2] Derivative: 20x³ + 12x² + 6x + 2
        assertTrue(Ex1.equals(Ex1.derivative(p), expected)); /// [3] Verify derivative works for large polynomials
    }





    //-------------------------end of tests derivative--------------------------------
	@Test
	/** 
	 * Tests the parsing of a polynom in a String like form.
	 */
	public void testFromString() {
		double[] p = {-1.1,2.3,3.1}; // 3.1X^2+ 2.3x -1.1
		String sp2 = "3.1x^2 +2.3x -1.1";
		String sp = Ex1.poly(p);
		double[] p1 = Ex1.getPolynomFromString(sp);
		double[] p2 = Ex1.getPolynomFromString(sp2);
		boolean isSame1 = Ex1.equals(p1, p);
		boolean isSame2 = Ex1.equals(p2, p);
		if(!isSame1) {fail();}
		if(!isSame2) {fail();}
		assertEquals(sp, Ex1.poly(p1));
	}
	@Test
	/**
	 * Tests the equality of pairs of arrays.
	 */
	public void testEquals() {
		double[][] d1 = {{0}, {1}, {1,2,0,0}};
		double[][] d2 = {Ex1.ZERO, {1+ Ex1.EPS/2}, {1,2}};
		double[][] xx = {{-2* Ex1.EPS}, {1+ Ex1.EPS*1.2}, {1,2, Ex1.EPS/2}};
		for(int i=0;i<d1.length;i=i+1) {
			assertTrue(Ex1.equals(d1[i], d2[i]));
		}
		for(int i=0;i<d1.length;i=i+1) {
			assertFalse(Ex1.equals(d1[i], xx[i]));
		}
	}

	@Test
	/**
	 * Tests is the sameValue function is symmetric.
	 */
	public void testSameValue2() {
		double x1=-4, x2=0;
		double rs1 = Ex1.sameValue(po1,po2, x1, x2, Ex1.EPS);
		double rs2 = Ex1.sameValue(po2,po1, x1, x2, Ex1.EPS);
		assertEquals(rs1,rs2, Ex1.EPS);
	}
	@Test
	/**
	 * Test the area function - it should be symmetric.
	 */
	public void testArea() {
		double x1=-4, x2=0;
		double a1 = Ex1.area(po1, po2, x1, x2, 100);
		double a2 = Ex1.area(po2, po1, x1, x2, 100);
		assertEquals(a1,a2, Ex1.EPS);
}
	@Test
	/**
	 * Test the area f1(x)=0, f2(x)=x;
	 */
	public void testArea2() {
		double[] po_a = Ex1.ZERO;
		double[] po_b = {0,1};
		double x1 = -1;
		double x2 = 2;
		double a1 = Ex1.area(po_a,po_b, x1, x2, 1);
		double a2 = Ex1.area(po_a,po_b, x1, x2, 2);
		double a3 = Ex1.area(po_a,po_b, x1, x2, 3);
		double a100 = Ex1.area(po_a,po_b, x1, x2, 100);
		double area =2.5;
		assertEquals(a1,area, Ex1.EPS);
		assertEquals(a2,area, Ex1.EPS);
		assertEquals(a3,area, Ex1.EPS);
		assertEquals(a100,area, Ex1.EPS);
	}
	@Test
	/**
	 * Test the area function.
	 */
	public void testArea3() {
		double[] po_a = {2,1,-0.7, -0.02,0.02};
		double[] po_b = {6, 0.1, -0.2};
		double x1 = Ex1.sameValue(po_a,po_b, -10,-5, Ex1.EPS);
		double a1 = Ex1.area(po_a,po_b, x1, 6, 8);
		double area = 58.5658;
		assertEquals(a1,area, Ex1.EPS);
	}
}
