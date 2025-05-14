import galois
import random
import numpy as np

class RScoderKaya(object):
    
    def __init__(self,n,k):
        """
        Initialize the Reed-Solomon encoder/decoder.

        This method initializes the Reed-Solomon encoder/decoder with a given codeword length and message length.

        Args:
            n (int): The length of the codeword.
            k (int): The length of the message.

        Raises:
        ValueError: If n or k are not positive, if n is not less than 256, or if k is not less than n.
        """
            
        GF=galois.GF(2**8, irreducible_poly="x^8 + x^4 + x^3 + x + 1",compile="jit-calculate")
        if n < 0 or k < 0:
            raise ValueError("n and k must be positive")
        if not n < 256:
            raise ValueError("n must be at most 255")
        if not k < n:
            raise ValueError("Codeword length n must be greater than message"
                    " length k")
        self.n = n
        self.k = k

        # Generate the generator polynomial for RS codes
        # g(x) = (x-α^1)(x-α^2)...(x-α^(n-k))
        # α is 3, a generator for GF(2^8)
        g=GF(1)
        for alpha in range(1,n-k+1):
            p = galois.Poly([1, GF(3)**alpha],field=GF)
            g = g * p
        self.g = g

        # h(x) = (x-α^(n-k+1))...(x-α^n)
        h = GF([1])
        for alpha in range(n-k+1,n+1):
            p = galois.Poly([1, GF(3)**alpha],field=GF)
            h = h * p
        self.h = h
        
 
    def base_convert(self,i, b):
        result = []
        while i > 0:
            result.insert(0, i % b)
            i = i // b
    
        while len(result)!=4:
            result.insert(0,0)
    
        return result

    def encode(self, message):
        """
            Encode a message using Reed-Solomon encoding.

            This function takes a message, encodes it using Reed-Solomon encoding, and returns the encoded message.

            Args:
                message (str): The message to encode.

            Returns:
                list: The encoded message as a list of coefficients.

            Raises:
                ValueError: If the length of the message is greater than k*4.
        """
        GF=galois.GF(2**8, irreducible_poly="x^8 + x^4 + x^3 + x + 1",compile="jit-calculate")
        n = self.n
        k = self.k
        if (len(message))>k*4:
            raise ValueError("Message length is max %d. Message was %d" % (k,
                len(message)))
        
        mtemp=0
        coefs=[]
        
        m = galois.Poly(message,field=GF)

        # Shift polynomial up by n-k by multiplying by x^(n-k)
        xshift=[1]
        for i in range(n-k): #MAYBE AN ERROR HERE POLYNOME ONLY GOES TO 254 WITHOUT THE 1
            xshift.append(0)
        mprime = m * galois.Poly(xshift,field=GF)

        # mprime = q*g + b for some q
        # so let's find b:
        b = divmod(mprime, self.g)[1]
        c = mprime - b
        ccof=c.coeffs
        coefstemp=[]
        for i in range(len(ccof)):
            coefstemp.append(int(ccof[i]))
            
        return coefstemp
    

    def verify(self,code):
        """
        Verify a Reed-Solomon code.

        This method takes a Reed-Solomon code and verifies it by checking if the remainder of the division of the code by the generator polynomial is zero.

        Args:
            code (list): The Reed-Solomon code to verify.

        Returns:
            bool: True if the code is valid, False otherwise.
        """
        #Maybe do it with h but haven't figured out how to quickly make an polynomial with max rank
        GF=galois.GF(2**8, irreducible_poly="x^8 + x^4 + x^3 + x + 1",compile="jit-calculate")
        g = self.g
        c = galois.Poly(code,field=GF)
        return c % g== GF(0)
    
    
    def syndromes(self,code):
        """
        Calculate the syndromes of a Reed-Solomon code.

        This method takes a Reed-Solomon code and calculates its syndromes by evaluating the code at different points in the finite field.

        Args:
            code (galois.Poly): The Reed-Solomon code to calculate the syndromes for.

        Returns:
            galois.Poly: The syndromes as a polynomial in the finite field.
        """
        GF=galois.GF(2**8, irreducible_poly="x^8 + x^4 + x^3 + x + 1",compile="jit-calculate")
        n = self.n
        k = self.k

        s=[GF(0)]
        for l in range(1, n-k+1):
            s.append( int(code(GF(3)**l)) )

        stemp=[]
        for lmi in range(len(s)-1):
            stemp.append(s[len(s)-1-lmi])
        stemp.append(int(GF(0)))
        sz = galois.Poly(stemp,field=GF) #Maybe not reversed?

        return sz


    def BMprosses(self,synd):
        """
        Perform the Berlekamp-Massey process for error correction.

        This method takes a syndromes polynomial and performs the Berlekamp-Massey process to find the error locator polynomial (sigma) and the error evaluator polynomial (omega).

        Args:
            synd (galois.Poly): The syndromes polynomial.

        Returns:
            tuple: The error locator polynomial (sigma) and the error evaluator polynomial (omega).
        """
        GF=galois.GF(2**8, irreducible_poly="x^8 + x^4 + x^3 + x + 1",compile="jit-calculate")
        n = self.n
        k = self.k

        # Initialize:
        sigma =  [galois.Poly([1],field=GF)]
        omega =  [galois.Poly([1],field=GF)]
        tao =    [galois.Poly([1],field=GF)]
        gamma =  [galois.Poly([0],field=GF)]
        D =      [0]
        B =      [0]

        # Polynomial constants:
        ONE = galois.Poly([1],field=GF)
        ZERO = galois.Poly([0],field=GF)
        Z = galois.Poly([1,0],field=GF)

        for l in range(0, n-k):
            # First find Delta, the non-zero coefficient of z^(l+1) in
            # (1 + s) * sigma[l]
            # This delta is valid for l (this iteration) only
            Deltatemp = ( (ONE + synd) * sigma[l] )

            delcoef=Deltatemp.coeffs
            delcoeftemp=[]
            for lmis in range(len(delcoef)):
                delcoeftemp.append(delcoef[len(delcoef)-1-lmis])
            Delta=GF(delcoeftemp[l+1])

            # Make it a polynomial of degree 0
            #Deltapol = galois.Poly([Delta],field=GF)

            sigma.append( sigma[l] - Delta * Z * tao[l] )
            omega.append( omega[l] - Delta * Z * gamma[l] )
            
            #if Delta==0:
            #    Derltainv=0
            #else:
            #    Derltainv=Delta**-1
            
            # Now compute the next tao and gamma
            # There are two ways to do this
            if Delta == ZERO or 2*D[l] > (l+1):
                # Rule A
                D.append( D[l] )
                B.append( B[l] )
                tao.append( Z * tao[l] )
                gamma.append( Z * gamma[l] )

            elif Delta != ZERO and 2*D[l] < (l+1):
                # Rule B
                D.append( l + 1 - D[l] )
                B.append( 1 - B[l] )
                tao.append( sigma[l] * Delta**-1)
                gamma.append( omega[l] * Delta**-1)
            elif 2*D[l] == (l+1):
                if B[l] == 0:
                    # Rule A (same as above)
                    D.append( D[l] )
                    B.append( B[l] )
                    tao.append( Z * tao[l] )
                    gamma.append( Z * gamma[l] )

                else:
                    # Rule B (same as above)
                    D.append( l + 1 - D[l] )
                    B.append( 1 - B[l] )
                    tao.append( sigma[l] * Delta**-1)
                    gamma.append( omega[l] * Delta**-1)
            else:
                raise Exception("Code shouldn't have gotten here")


        return sigma[-1], omega[-1]
    

    def Csearch(self, sigma):
        """
        Search for error locations in a Reed-Solomon code.

        This method takes an error locator polynomial (sigma) and finds the locations of errors in the Reed-Solomon code by finding the roots of the polynomial.

        Args:
            sigma (galois.Poly): The error locator polynomial.

        Returns:
            tuple: The error locations (X) and their corresponding indices (j).
        """
        GF=galois.GF(2**8, irreducible_poly="x^8 + x^4 + x^3 + x + 1",compile="jit-calculate")
        X = []
        j = []
        p = GF([3])
        sigzeros=sigma.roots()
        for l in range(1,256):
            # These evaluations could be more efficient, but oh well
            #x=p**l
            #m=len(sigzeros)
            #while m != 0:
            if sigma(GF(3)**l) == 0:
                X.append( GF(3)**(-l) )
                # This is different than the notes, I think the notes were in error
                # Notes said j values were just l, when it's actually 255-l
                j.append(255 - l)

                #m=m-1 
        return X, j
    
    def Fformular(self, omega, X):
        """
        Compute the error values in a Reed-Solomon code.

        This method takes an error evaluator polynomial (omega) and the error locations (X), and computes the error values (Y) using Forney's formula.

        Args:
            omega (galois.Poly): The error evaluator polynomial.
            X (list): The error locations.

        Returns:
            list: The error values (Y).
        """
        # XXX Is floor division okay here? Should this be ceiling?
        
        GF=galois.GF(2**8, irreducible_poly="x^8 + x^4 + x^3 + x + 1",compile="jit-calculate")
        s = (self.n - self.k) // 2

        Y = []

        for l, Xl in enumerate(X):
            # Compute the first part of Yl
            Yl = Xl**s
            Yl *= omega( Xl** -1 )
            Yl *= Xl** -1

            # Compute the sequence product and multiply its inverse in
            prod = GF(1)
            for ji in range(s):
                if ji == l:
                    continue
                if ji < len(X):
                    Xj = X[ji]
                else:
                    Xj = GF(0)
                prod = prod * (Xl - Xj)
            if prod==0:
                Yl = 0
            else:
                Yl = Yl * prod** -1

            Y.append(Yl)
        return Y



    def decode(self, codegiv):
        """
        Decode a Reed-Solomon code.

        This method takes a Reed-Solomon code and decodes it. If the code is valid, it removes the parity bits and returns the decoded message. If the code is not valid, it performs error correction using the Berlekamp-Massey process and Forney's formula, and then returns the corrected message.

        Args:
            codegiv (list): The Reed-Solomon code to be decoded.

        Returns:
            list: The decoded message.
        """
        GF=galois.GF(2**8, irreducible_poly="x^8 + x^4 + x^3 + x + 1",compile="jit-calculate")
        n = self.n
        k = self.k

        if self.verify(codegiv):
            # The last n-k bytes are parity
            lastbits=n-k
            final=[]
            for j in range(len(codegiv)):
                final.append(codegiv[j])
            for i in range(lastbits):
                final.pop()
            len(final)
            quartfin=[]
            for i in range(len(final)):
               quartfin.extend(self.base_convert(int(final[i]),4))
            return quartfin


        c=galois.Poly(codegiv,field=GF)
        sz = self.syndromes(c)
        sigma, omega = self.BMprosses(sz)
        X, j = self.Csearch(sigma)
        Y = self.Fformular(omega, X)

        # Put the error and locations together to form the error polynomial
        Elist = [0]
        for i in range(255):
            if i in j:
                Elist.append(Y[j.index(i)])
            else:
                Elist.append(GF(0))
        Elistrev=[]
        for lmi in range(len(Elist)-1):
            Elistrev.append(Elist[len(Elist)-1-lmi])
        E = galois.Poly(Elistrev,field=GF)
        final = c - E

        finalcoefs=final.coeffs
        finalfcoefs=[]
        for j in range(len(finalcoefs)):
            finalfcoefs.append(finalcoefs[j])
        for i in range(n-k):
            finalfcoefs.pop()
        coefsf=[]
        for i in range(len(finalfcoefs)):
           coefsf.extend(self.base_convert(int(finalfcoefs[i]),4))

        return coefsf





