from math import comb
def sign_test_power(n, p, alpha=0.05):
    # one-sided exact sign test on n discordant pairs; P(success)=p under the alternative
    # critical k: smallest k with P(X>=k | 0.5) <= alpha
    k = next(k for k in range(n+1) if sum(comb(n,j) for j in range(k,n+1))/2**n <= alpha)
    return k, sum(comb(n,j)*p**j*(1-p)**(n-j) for j in range(k,n+1))
print('n_discordant  crit_k  power@p=0.65  p=0.70  p=0.75  p=0.80  p=0.90')
for n in (5,8,10,12,15,20,25,30,40,50,60,80):
    k,_=sign_test_power(n,0.5)
    print(f'{n:5d} {k:9d}', '  '.join(f'{sign_test_power(n,p)[1]:.2f}' for p in (0.65,0.70,0.75,0.80,0.90)))
