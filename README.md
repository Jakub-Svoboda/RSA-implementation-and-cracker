# RSA-implementation-and-cracker
The implementation of a RSA encryption protocol. Also includes a naive cracker for RSA. Implemented for KRY class at FIT VUT.

For key generation, the p and q prime numbers are found using the solovay-strassen algorithm. 

For factorization, the Extended Eucliedean method is used.

The cracking method for RSA is implemented using naive factorization algorithm, where the algorithm searches for a number
which when used as a denominator on N results in zero. For this, the algorithm simple tries every number larger than 3 and
increments by 2 after each unsuccessfull try.
