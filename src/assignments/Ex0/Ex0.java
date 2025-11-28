package assignments.Ex0;
/**
 * This class is a basis for Ex0 (your first assigment),
 * The definition of the Ex0 can be found here: https://docs.google.com/document/d/1UtngN203ttQKf5ackCnXs4UnbAROZWHr/edit?usp=sharing&ouid=113711744349547563645&rtpof=true&sd=true
 * You are asked to complete the functions below and may add additional functions if needed.

 */
public class Ex0 {
    public final static long ID = 330920141;  // Do update your ID here

    /// סעיף 1


    /**
     *
     * isPrime(n)
     *    if( n!=2 && n % 2 == 0 || n!=3 && n % 3 == 0 || n<2 )
     *            return F
     *      for long i = 5; Math.sqrt <= n; i += 6) {
     *                  if (n % i == 0 OR n % (i + 2) == 0)
     *                      return F
     *              return T
     */


    public static boolean isPrime(long n) {
        /// if not 2, not even
        /// or not 3 and divided by 3
        /// or if n is negative num or 0 or 1(that smaller than 2)
        if (n!=2 && n % 2 == 0 || n!=3 && n % 3 == 0 || n<2)
            return false;
        for (long i = 5; i*i <= n; i += 6) {
            if (n % i == 0 || n % (i + 2) == 0)
                return false;
        }
        return true;
    }


    /// סעיף 2


    /**
     * f(int num)
     * int prime=num+1
     * bool flag=F
     * while(flag==F)
     *      if(isprime(prime))
     *          flag=T
     *      else
     *          prime++
     *return prime
     */

    /// FNP= first next prime( that is bigger than s)
    public static long FNP(long num) {
        long prime = num + 1;
        /// if not prime, then we will look for the next close parameter if it is a prime number
        while (!isPrime(prime)) {
            prime++;
        }
        return prime;
    }

    /**
     * getprimepair(s,n)
     * prime = FNP(s)
     * ans=-1
     * while (T)
     *      if (isPrime(prime + n))
     *                 break
     *      else
     *                 prime = FNP(prime)
     *return prime
     */

    /// a function that returns prime num only if the next num p2=p1+n is also prime

    public static long getPrimePair(long s, long n) {
        long prime = FNP(s);
        long ans=-1;
        while (true) {
            ///checks if p2=p1+n
            if (isPrime(prime + n)) {
                /// breaks from the while
                break;
            }
            /// if not, then we will look for the next prime num that will answer our condition
            else
                prime = FNP(prime);

        }
        return prime;
    }


    /// סעיף 3


    /**
     * getclosestprimepair(s,n)
     *         ans = -1;
     *         if (n >= 2 AND n % 2 == 0)
     *             i = s
     *             if (i % 2 == 0) i+1
     *             while (T)
     *                 if (isPrime(i + n))
     *                     if (isPrime(i))
     *                         flag = F
     *                         for (long j = i + 2; j < n + i; j += 2)
     *                             if (isPrime(j))
     *                                 flag = T
     *                                 break

     *                         if (!flag)
     *                             return i

     *                     i = i + n
     *                  else
     *                     i += 2
     *         return ans
     */


    public static long getClosestPrimePair(long s, long n) {
        /// Default return value if no pair is found
        long ans = -1;

        /// Check if n is even and bigger or equal to 2
        if (n >= 2 && n % 2 == 0) {
            /// Start from s
            long i = s;

            /// Ensure i is odd
            if (i % 2 == 0) i++;

            /// Loop until a valid prime pair is found
            while (true) {
                /// Check if (i + n) is prime
                if (isPrime(i + n)) {
                    /// Check if i is prime
                    if (isPrime(i)) {
                        /// Flag to check if another prime exists between i and i+n
                        boolean flag = false;

                        for (long j = i + 2; j < n + i; j += 2) {
                            /// If another prime is found, flag=T and break
                            if (isPrime(j)) {
                                flag = true;
                                break;
                            }
                        }

                        /// If no other prime exists
                        if (!flag) {
                            return i;
                        }
                    }

                    /// Move i forward by n
                    i = i + n;
                } else {
                    /// If (i + n) is not prime, move i by 2
                    i += 2;
                }
            }
        }

        /// Return -1 if no pair found
        return ans;
    }



    /// סעיף 4
    /**
     * getMthclosestprimepair(m,n)
     * s=2
     * count=0
     * result=0
     * while(count<m)
     *      candidate = getClosestPrimePair(s, n)
     *      result = candidate
     *      count++;
     *      s = candidate + n + 1
     *return result
     */


    public static long getMthClosestPrimePair(int m, long n) {
        int count = 0;
        long s = 2;
        long result = 0;
        /// runs untill it finds the pair in the order
        while (count <= m) {
            long candidate = getClosestPrimePair(s, n);
            result = candidate;
            count++;
            /// continue after the pair that we found
            s = candidate + 1;
        }
        return result;
    }

}








