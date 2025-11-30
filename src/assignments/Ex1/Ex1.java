package assignments.Ex1;

/**
 * Introduction to Computer Science 2026, Ariel University,
 * Ex1: arrays, static functions and JUnit
 * https://docs.google.com/document/d/1GcNQht9rsVVSt153Y8pFPqXJVju56CY4/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 *
 * This class represents a set of static methods on a polynomial functions - represented as an array of doubles.
 * The array {0.1, 0, -3, 0.2} represents the following polynomial function: 0.2x^3-3x^2+0.1
 * This is the main Class you should implement (see "add your code below")
 *
 * @author boaz.benmoshe

 */
public class Ex1 {
    /**
     * Epsilon value for numerical computation, it serves as a "close enough" threshold.
     */
    public static final double EPS = 0.001; // the epsilon to be used for the root approximation.
    /**
     * The zero polynomial function is represented as an array with a single (0) entry.
     */
    public static final double[] ZERO = {0};

    /**
     * Computes the f(x) value of the polynomial function at x.
     *
     * @param poly - polynomial function
     * @param x
     * @return f(x) - the polynomial function value at x.
     */
    public static double f(double[] poly, double x) {
        double ans = 0;
        for (int i = 0; i < poly.length; i++) {
            double c = Math.pow(x, i);
            ans += c * poly[i];
        }
        return ans;
    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an epsilon eps.
     * This function computes an x value (x1<=x<=x2) for which |p(x)| < eps,
     * assuming p(x1)*p(x2) <= 0.
     * This function should be implemented recursively.
     * @param p   - the polynomial function
     * @param x1  - minimal value of the range
     * @param x2  - maximal value of the range
     * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
     * @return an x value (x1<=x<=x2) for which |p(x)| < eps.
     */
    public static double root_rec(double[] p, double x1, double x2, double eps) {
        double f1 = f(p, x1);
        double x12 = (x1 + x2) / 2;
        double f12 = f(p, x12);
        if (Math.abs(f12) < eps) {
            return x12;
        }
        if (f12 * f1 <= 0) {
            return root_rec(p, x1, x12, eps);
        } else {
            return root_rec(p, x12, x2, eps);
        }
    }

    /**
     * This function computes a polynomial representation from a set of 2D points on the polynom.
     * The solution is based on: //	http://stackoverflow.com/questions/717762/how-to-calculate-the-vertex-of-a-parabola-given-three-points
     * Note: this function only works for a set of points containing up to 3 points, else returns null.
     *
     * @param xx
     * @param yy
     * @return an array of doubles representing the coefficients of the polynom.
     */
    public static double[] PolynomFromPoints(double[] xx, double[] yy) {
        double[] ans = null;
        int lx = xx.length;
        int ly = yy.length;
        if (xx != null && yy != null && lx == ly && lx > 1 && lx < 4) {
            /** add you code below

             /////////////////// */
        }
        return ans;
    }

    /**
     * Two polynomials functions are equal if and only if they have the same values f(x) for n+1 values of x,
     * where n is the max degree (over p1, p2) - up to an epsilon (aka EPS) value.
     *
     * @param p1 first polynomial function
     * @param p2 second polynomial function
     * @return true iff p1 represents the same polynomial function as p2.
     */
    public static boolean equals(double[] p1, double[] p2) {
        boolean ans = true;
        /** add you code below

         /////////////////// */
        return ans;
    }

    /**
     * Computes a String representing the polynomial function.
     * For example the array {2,0,3.1,-1.2} will be presented as the following String  "-1.2x^3 +3.1x^2 +2.0"
     *
     * @param poly the polynomial function represented as an array of doubles
     * @return String representing the polynomial function:
     */
    public static String poly(double[] poly) {
        String ans = "";
        if (poly.length == 0) {
            ans = "0";
        } else {
            /** add you code below

             /////////////////// */
        }
        return ans;
    }

    /**
     * Given two polynomial functions (p1,p2), a range [x1,x2] and an epsilon eps. This function computes an x value (x1<=x<=x2)
     * for which |p1(x) -p2(x)| < eps, assuming (p1(x1)-p2(x1)) * (p1(x2)-p2(x2)) <= 0.
     *
     * @param p1  - first polynomial function
     * @param p2  - second polynomial function
     * @param x1  - minimal value of the range
     * @param x2  - maximal value of the range
     * @param eps - epsilon (positive small value (often 10^-3, or 10^-6).
     * @return an x value (x1<=x<=x2) for which |p1(x) - p2(x)| < eps.
     */
    public static double sameValue(double[] p1, double[] p2, double x1, double x2, double eps) {
        double ans = x1;
        /** add you code below

         /////////////////// */
        return ans;
    }

    /**
     * Given a polynomial function (p), a range [x1,x2] and an integer with the number (n) of sample points.
     * This function computes an approximation of the length of the function between f(x1) and f(x2)
     * using n inner sample points and computing the segment-path between them.
     * assuming x1 < x2.
     * This function should be implemented iteratively (none recursive).
     *
     * @param p                - the polynomial function
     * @param x1               - minimal value of the range
     * @param x2               - maximal value of the range
     * @param numberOfSegments - (A positive integer value (1,2,...).
     * @return the length approximation of the function between f(x1) and f(x2).
     */
    public static double length(double[] p, double x1, double x2, int numberOfSegments) {
        double ans = x1;
        /** add you code below

         /////////////////// */
        return ans;
    }

    /**
     * Given two polynomial functions (p1,p2), a range [x1,x2] and an integer representing the number of Trapezoids between the functions (number of samples in on each polynom).
     * This function computes an approximation of the area between the polynomial functions within the x-range.
     * The area is computed using Riemann's like integral (https://en.wikipedia.org/wiki/Riemann_integral)
     *
     * @param p1                - first polynomial function
     * @param p2                - second polynomial function
     * @param x1                - minimal value of the range
     * @param x2                - maximal value of the range
     * @param numberOfTrapezoid - a natural number representing the number of Trapezoids between x1 and x2.
     * @return the approximated area between the two polynomial functions within the [x1,x2] range.
     */
    public static double area(double[] p1, double[] p2, double x1, double x2, int numberOfTrapezoid) {
        double ans = 0;
        /** add you code below

         /////////////////// */
        return ans;
    }

    /**
     * This function computes the array representation of a polynomial function from a String
     * representation. Note:given a polynomial function represented as a double array,
     * getPolynomFromString(poly(p)) should return an array equals to p.
     *
     * @param p - a String representing polynomial function.
     * @return
     */
    public static double[] getPolynomFromString(String p) {
        double[] ans = ZERO;//  -1.0x^2 +3.0x +2.0
        /** add you code below

         /////////////////// */
        return ans;
    }

    /**
     * This function computes the polynomial function which is the sum of two polynomial functions (p1,p2)

     *public static void replace(double[], double[]) /// the function is void because it (will\will not) just swap the references of the two polynoms
     *    /// if (p1.length > p2.length || p1.length == p2.length) then continue as it is
     *    else if(p2.length > p1.length)
     *      temp = p1
     *      p1 = p2
     *      p2 = temp
     *
     * add(p1[], p2[])
     *     if (p1 == null && p2 == null)
     *     return ZERO /// return {0}
     *     if (p1 == null && p2!=null )
     *          ans = p2
     *     if (p2 == null && p1!=null )
     *          ans = p1
     *      len=p1.length-p2.length
     *      ans = new double[p1.length]
     *      for(i=0; i<len; i++)
     *          ans[i] = p1[i]
     *      for(i=0; i<p2.length; i++)
     *          ans[i+len] = p2[i] + p1[i+len]
     *      return ans
     * @return
     */
    public static double[] compact(double[] p) {
        int i=0;
        while(i<=p.length-1 && p[i]!=0 )
        {
            i++;
        }
        double[] ans = new double[i];
        for (int j=0; j<ans.length; j++)
        {
            ans[j] = p[j];
        }
        return ans;
    }
    public static double[] add(double[] p1, double[] p2) {

       /// first we will compact both polynoms to make sure that there are no unnecessary zeros at the beginning of the polynoms
       if (p1 != null && p1.length>0)
            p1=compact(p1);
        if (p2 != null && p2.length>0)
            p2=compact(p2);
        int len1, len2;
        /// find the lengths of both polynoms
        if(p1 != null && p1.length>0)
            len1= p1.length;
        else {
            len1=0;
        }
        if (p2 != null && p2.length>0)
            len2= p2.length;
        else {
            len2=0;
        }
        /// to make sure that p1 is the longer polynom
        /// if (p1.length > p2.length || p1.length == p2.length) then continue as it is, else swap p1 and p2, so p1 will always be the longer polynom
        if (len2 > len1) {
            double[] temp = p1;
            p1 = p2;
            p2 = temp;
            /// swap lengths as well
            int tempLen = len1;
            len1 = len2;
            len2 = tempLen;
        }
        double[] ans;
        /// if both polynoms are empty return ZERO
        if ((p1 == null && len1==0) && (p2 == null || len2==0))
          return ZERO; /// return {0}
        if ((p1 == null || len1>0) && (p2!=null && len2>0)) {
            ans = p2; /// if p1 is empty return p2 bcause p2 + 0 = p2
            return ans;
        }
        if ((p2 == null || len2==0) && (p1!=null && len1>0)) {
            ans = p1; /// if p2 is empty return p1 because p1 + 0 = p1
        return ans;
        }
        int len=len1-len2; /// here we find the difference between the lengths of p1 and p
        ans = new double[len1]; /// because p1 is longer so it has the max elements that ans will have when we add p1+p2
        if(len!=0) { /// check if len is not zero because if it is zero that means that both polynoms have the same length so there is no need to add the first len elements from p1 to ans
            for (int i = 0; i < len; i++) { /// we will run this loop len times because we want to add all the elements from the first element of p1 which dont have a the same power as in p2
                ans[i] = p1[i]; /// put into ans the first len elements from p1 because p2 has no elements in these pow
            }
        }
        for(int i=0; i<len2; i++) { ///we will run this loop p2.length times because we want to add all the elements from the first element of p2
            ans[i + len] = p2[i] + p1[i + len]; /// add the rest of the elements from p1 and p2
        }

        return ans;
    }

    /**
     * This function computes the polynomial function which is the multiplication of two polynoms (p1,p2)
     mul[](p1[],p2[]) {

     if (p1 != null && p1.length>0)
        p1=compact(p1)
     if (p2 != null && p2.length>0)
     p2=compact(p2)
        len1, len2;

     if(p1 != null && p1.length>0)
         len1= p1.length
         else
         len1 = 0

     if (p2 != null && p2.length>0)
         len2= p2.length
         else
         len2 = 0

     if (len2 > len1)
         temp[] = p1
         p1 = p2
         p2 = temp
         tempLen = len1
         len1 = len2
         len2 = tempLen


     int len=len1-len2
     ans[] = new double[len1]

     if(p1==null || p2==null || len1 == 0 || len2 == 0)
        ans=null
        return ans


     if(len!=0)
         int index1=len1-1
         for (int i=len-1; i >= 0; i--)
         ans[index1] = p1[index1]
         index1--
         i=len-1
         j=len2-1
         index=len2-1
         while(i>=0 && j>=0)
         ans[index]=p1[i]*p2[j]
         i--
         j--
         index--


     else
     if (len1==1 && p1[0] ==1)
        ans= p2;

     else if (len2==1 && p2[0] ==1)
        ans = p1

     else if(p1==ZERO || p2==ZERO)
        ans= 0

     return ans
     */

}

    public static double[] mul(double[] p1, double[] p2) {

        /// first we will compact both polynoms to make sure that there are no unnecessary zeros at the beginning of the polynoms
        if (p1 != null && p1.length>0)
            p1=compact(p1);
        if (p2 != null && p2.length>0)
            p2=compact(p2);
        int len1, len2;
        /// find the lengths of both polynoms
         if(p1 != null && p1.length>0)
             len1= p1.length;
         else {
             len1 = 0;
         }
         if (p2 != null && p2.length>0)
             len2= p2.length;
         else {
             len2 = 0;
         }
            /// to make sure that p1 is the longer polynom
            /// if (p1.length > p2.length || p1.length == p2.length) then continue as it is, else swap p1 and p2, so p1 will always be the longer polynom
        if (len2 > len1) {
                double[] temp = p1;
                p1 = p2;
                p2 = temp;
                /// swap lengths as well
                int tempLen = len1;
                len1 = len2;
                len2 = tempLen;
        }
        ///double[] ans;
        int len=len1-len2; /// here we find the difference between the lengths of p1 and p
        double[] ans = new double[len1]; /// because p1 is longer so it has the max elements that ans will have when we add p1+p2

             if(p1==null || p2==null || len1 == 0 || len2 == 0)   /// if one of the polynoms is null then the answer is null because multiplying by null is not defined
             {
                ans=null;
                return ans;
             }

             if(len!=0) {  /// check if len is not zero because if it is zero that means that both polynoms have the same length so there is no need to add the first len elements from p1 to ans
               int index1=len1-1;
               for (int i=len-1; i >= 0; i--) {  /// we will run this loop len times because we want to mul all the elements from the first element of p1 which dont have a the same power as in p2
                    ans[index1] = p1[index1];  /// put into ans the first len elements from p1 because p2 has no elements in these pow
                    index1--;
               }
               int i=len-1;
               int j=len2-1;
               int index=len2-1;
               while(i>=0 && j>=0) { /// we will run this loop p2.length times because we want to mul all the elements from the first element of p2
                   ans[index]=p1[i]*p2[j];
                    i--;
                    j--;
                    index--;
               }

           }
           else  /// if both polynoms have the same length and one of them is {1} then the answer is the other polynom
           {
               if (len1==1 && p1[0] ==1) {   /// if p1 is {1} then the answer is p2 because 1*p2 = p2
                   ans= p2;
               }
               else if (len2==1 && p2[0] ==1) {  /// if p2 is {1} then the answer is p1 because 1*p1 = p1
                   ans = p1;
               }
               else if(p1==ZERO || p2==ZERO) {   /// if one of the polynoms is ZERO then the answer is ZERO because 0*p1 =0
                   ans= ZERO;
               }
           }
        return ans;
    }









    /**
     * This function computes the derivative of the p0 polynomial function.
     *
     * derivative(po[] /// ={1,2,3} represents 3x^2 +2x +1)
     *       ans = ZERO /// ans = {0}
     *       if (po != null && po.length > 1) /// T && T ==T
     *             len = po.length /// len = 3
     *             ans = new double[len - 1] /// ans = new double[2] /// ans = {0,0}
     *             for (i = 0; i < ans.length; i=i+1)
     *                  ans[i] = po[i + 1] * (i + 1) /// ans[0] = po[1]*1 = 2*1 =2 , ans[1] = po[2]*2 = 3*2=6
     *       return ans /// ans = {2,6} represents 6x +2
     */
    public static double[] derivative(double[] po) {
        double[] ans = ZERO;
        /// check if po is not empty and has more than one element in it
        if (po != null && po.length > 1) {
            int len = po.length;
            /// create the answer array with length less by one because the derivative degree is less by one than the original polynom
            ans = new double[len - 1];
            /// compute the derivative using the power rule of derivatives ([a*x^n] = n*a*x^(n-1) )
            for (int i = 0; i < ans.length; i += 1) {
                /// puts into array ans the num from the place po[1] to the end multiplied by its degree (i+1)
                ans[i] = po[i + 1] * (i + 1);
            }
        }
        /// return the derivative polynom in a new array
        return ans;
    }
}
