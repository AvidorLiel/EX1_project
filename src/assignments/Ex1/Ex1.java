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
     *
     *
     * mid = (x1 + x2) / 2
     *  if |f(mid)| < eps:
     *      return mid
     *  if f(x1)*f(mid) <= 0:
     *      recurse on [x1, mid]
     *  else:
     *      recurse on [mid, x2]
     */
    public static double root_rec(double[] p, double x1, double x2, double eps) {
        double f1 = f(p,x1); /// compute f(x1)
        double x12 = (x1+x2)/2; /// compute the mid point
        double f12 = f(p,x12); /// compute f(mid)
        if (Math.abs(f12)<eps) /// check if |f(mid)| < eps
        {
            return x12; /// return mid
        }
        if(f12*f1<=0) /// check if f(x1)*f(mid) <= 0
        {
            return root_rec(p, x1, x12, eps); /// recurse on [x1, mid]
        }
        else
        {
            return root_rec(p, x12, x2, eps); /// recurse on [mid, x2]
        }
    }



    /**
     * This function computes the array representation of a polynomial function from a String
     * representation. Note:given a polynomial function represented as a double array,
     * getPolynomFromString(poly(p)) should return an array equals to p.
     *
     * @param p - a String representing polynomial function.
     * @return
     *
     *  split string into terms by +/-
     *  parse coefficient and power for each term
     *  accumulate into array
     */
    public static double[] getPolynomFromString(String p) {
        if (p == null || p.trim().isEmpty()) /// check for null or empty string
            return new double[]{0};

        p = p.replace(" ", ""); /// remove spaces

        java.util.List<String> terms = new java.util.ArrayList<>();
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < p.length(); i++) { /// split string into terms by +/-
            char c = p.charAt(i);
            if ((c == '+' || c == '-') && sb.length() > 0) { /// new term
                terms.add(sb.toString());
                sb.setLength(0);
            }
            sb.append(c);
        }
        if (sb.length() > 0) terms.add(sb.toString()); /// add last term

        java.util.Map<Integer, Double> coeffMap = new java.util.HashMap<>();
        int maxPower = 0;

        for (String term : terms) { /// parse coefficient and power for each term
            if (term.isEmpty()) continue;

            double sign = 1.0;
            if (term.charAt(0) == '-') { /// handle sign
                sign = -1.0;
                term = term.substring(1);
            } else if (term.charAt(0) == '+') { /// handle sign
                term = term.substring(1);
            }

            double coeff = 0.0;
            int power = 0;

            if (term.contains("x")) { /// parse x term
                String[] parts = term.split("x");
                if (parts[0].isEmpty()) /// handle coefficient
                    coeff = 1.0;
                else /// handle coefficient
                    coeff = Double.parseDouble(parts[0]);

                if (parts.length > 1 && parts[1].startsWith("^")) { /// handle power
                    power = Integer.parseInt(parts[1].substring(1));
                }
                else { /// handle power
                    power = 1;
                }
            } else { /// constant term
                coeff = Double.parseDouble(term);
                power = 0;
            }

            coeff *= sign; /// apply sign
            coeffMap.put(power, coeffMap.getOrDefault(power, 0.0) + coeff); /// accumulate into array
            maxPower = Math.max(maxPower, power); /// update max power
        }

        double[] ans = new double[maxPower + 1]; /// create result array
        for (int i = 0; i <= maxPower; i++) { /// fill result array
            ans[i] = coeffMap.getOrDefault(i, 0.0); /// get coefficient
        }
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
