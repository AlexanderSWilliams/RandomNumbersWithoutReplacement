
# coding: utf-8

# In[1]:

#Expressions and matricies that will be discussed
import numpy as np
from sympy import *
from fractions import *
def stateTElements(n):
    def f(i,j):
        if i==j:
            return j/n if type(j) == Symbol else Fraction(j,n)
        elif i-j==1:
            return 1-j/n if type(j) == Symbol else 1-Fraction(j,n)
        return 0
    return f
def stateDElement(n, M = 1):
    def f(i,j):
        if(i==j):
            return (j/n)**M if type(j) == Symbol else Fraction(j,n)**M
        return 0
    return f
def statePElement(n):
    def f(i,j):
        return (-1)**(i+j)*binomial(n,i)*binomial(i,j)
    return f
def statePIElement(n):
    def f(i,j):
            return binomial(i,j)/binomial(n,j) if type(i) == Symbol         else Fraction(int(binomial(i,j)), int(binomial(n,j)))
    return f
def stateT(n):
    return Matrix(n+1,n+1, [ stateTElements(n)(i,j) for i in range(n+1) for j in range(n+1)])
def stateD(n):
    return Matrix(n+1,n+1, [ stateDElement(n)(i,j) for i in range(n+1) for j in range(n+1)])
def stateP(n):
    return Matrix(n+1,n+1, [ statePElement(n)(i,j) for i in range(n+1) for j in range(n+1)])
def statePI(n):
    return Matrix(n+1,n+1, [ statePIElement(n)(i,j) for i in range(n+1) for j in range(n+1)])


# In[2]:

#What is stateT and why is stateT a transition matrix that represents choosing distinct elements without replacement?
stateT(5)


# In[3]:

#How to represent the probabilites of given states of M choices?
print(stateT(5)*Matrix(5+1,1,[1,0,0,0,0,0]))
print(stateT(5)**2*Matrix(5+1,1,[1,0,0,0,0,0]))
print(stateT(5)**3*Matrix(5+1,1,[1,0,0,0,0,0]))


# In[4]:

#What is the expectation of the number of distinct elements of M choice attempts?
result = Matrix(1,1+5,[0,1,2,3,4,5])*stateT(5)**3*Matrix(5+1,1,[1,0,0,0,0,0]);
print(result)
print(result.evalf(10))


# In[5]:

#How is expectation of the number of distinct elements significant to the algorithm?
print(Matrix(1,5+1,range(5+1))*stateT(5)**5*Matrix(5+1,1,[1 if i == 0 else 0 for i in range(5+1)]).evalf(10))
print(Matrix(1,10+1,range(10+1))*stateT(10)**10*Matrix(10+1,1,[1 if i == 0 else 0 for i in range(10+1)]).evalf(10))
print(Matrix(1,15+1,range(15+1))*stateT(15)**15*Matrix(15+1,1,[1 if i == 0 else 0 for i in range(15+1)]).evalf(10))
print(Matrix(1,20+1,range(20+1))*stateT(20)**20*Matrix(20+1,1,[1 if i == 0 else 0 for i in range(20+1)]).evalf(10))


# In[6]:

#What is the Jordan Canonical Form of stateT?
pprint(stateP(5)*stateD(5)*statePI(5))
pprint(stateD(5))
pprint(stateP(5)*statePI(5))


# In[7]:

#For n>1 stateP(n)*statePI(n) = I
n, i, j, k = symbols('n i j k', integer=True)
statePIElement(n)(i,k)*statePElement(n)(k,j) # need to sum from 0 <= k <= n


# In[8]:

#Anti-difference for stateP(n)*statePI(n)
antiDiffQIQ = (-1)**(j+k)*(i-k)*binomial(i,k)*binomial(k,j)/(i-j)
print(antiDiffQIQ - antiDiffQIQ.subs({k:k-1}) - statePIElement(n)(i,k)*statePElement(n)(k,j))
combsimp(antiDiffQIQ - antiDiffQIQ.subs({k:k-1}) - statePIElement(n)(i,k)*statePElement(n)(k,j))


# In[9]:

#Expression for the cells of stateP(n)*statePI(n) when i is not equal to j
antiDiffQIQ.subs({k:n}) - antiDiffQIQ.subs({k:0})


# In[10]:

#Expression for the cells of stateP(n)*statePI(n) when i is 
summation(statePIElement(n)(i,k)*statePElement(n)(k,i), (k, i, i)) 


# In[11]:

#What is the expression for the terms of stateP*stateD?
summation(statePElement(n)(i,k)*stateDElement(n)(j,j), (k,j,j))


# In[12]:

#What is the expression for the terms of stateP*stateD*statePI?
jordanTerms = summation(statePElement(n)(i,k)*stateDElement(n)(j,j), (k,j,j)).subs({j:k})*statePIElement(n)(k,j) 
jordanTerms # need to sum from 0 <= k <= n


# In[13]:

#What is the anti-difference for the terms of stateP*stateD*statePI?
antiDiffPDPI = ((-1)**(i+k)*((i-j)*k-j)*(i-k)*binomial(i,k)*binomial(k,j)/((i-j-1)*(i-j)))*binomial(n, i)/(n*binomial(n, j))
print(antiDiffPDPI - antiDiffPDPI.subs({k:k-1}) - jordanTerms)
combsimp(antiDiffPDPI - antiDiffPDPI.subs({k:k-1}) - jordanTerms)


# In[14]:

#Proof that the terms of stateP*stateD*statePI are zero when i != j and i != j + 1
combsimp(antiDiffPDPI.subs({k:n}) - antiDiffPDPI.subs({k:0}))


# In[15]:

#Expression for the i=j terms of stateP*stateD*statePI
iEqualsJTerms = jordanTerms.subs({i:j})
print(iEqualsJTerms)
iEqualsJTerms.subs({k:j})


# In[16]:

#Expression for the i=j+1 terms of stateP*stateD*statePI
iEqualsJPlusOneTerms = jordanTerms.subs({i:j+1})
print(iEqualsJPlusOneTerms)
combsimp(iEqualsJPlusOneTerms.subs({k:j}) + iEqualsJPlusOneTerms.subs({k:j+1}))


# In[17]:

#Expression for the i=j terms of stateP*stateD*statePI
A = symbols('A', integer=True)
jordanPowerTerms = summation(statePElement(n)(i,k)*stateDElement(n,n)(j,j), (k,j,j)).subs({j:k})*statePIElement(n)(k,j)
print(jordanPowerTerms)
summation(A*summation(jordanPowerTerms,(k,0,n)).subs({i:A,j:0}), (A,0,n))


# In[18]:

#Anti-difference for the expectation of the number of distinct elements.
antiDiffExpectation = ((-1)**(A+k)*((k-n)*A+k)*(A-n)*binomial(A,k)*binomial(n,A)/((k-n+1)*(k-n)))*(k/n)**n
print(combsimp(antiDiffExpectation - antiDiffExpectation.subs({A:A-1}) - A*jordanPowerTerms.subs({i:A,j:0})))
print(combsimp(antiDiffExpectation.subs({A:n}) - antiDiffExpectation.subs({A:0})))


# In[19]:

#Expression for the expectation of the number of distinct elements
print(summation(A*summation(jordanPowerTerms,(k,n-1,n)).subs({i:A,j:0}), (A,0,n)))
combsimp(summation(A*summation(jordanPowerTerms,(k,n-1,n)).subs({i:A,j:0}), (A,n-1,n)))


# In[20]:

#Sanity check of the expression
print((-n*((1 - 1/n)**n - 1)).subs({n:2}))
print((-n*((1 - 1/n)**n - 1)).subs({n:3}))
print(1+Fraction(1,2))
print(1+Fraction(2,3)+Fraction(2,3)*Fraction(1,3)+Fraction(1,3)*Fraction(2,3))


# In[21]:

#The constant of interest and its approximation.
print(simplify(limit((-n*((1 - 1/n)**n - 1))/n,n,oo)-(1-exp(-1))))
print((1-exp(-1)).evalf(10))


# In[22]:

#Empirical observation of the low variance of the expectation
from random import randrange
totalCount = 0
def getDistinctRandomNumbers(inclusiveMin, exclusiveMax, amount):
    global totalCount
    result = set()
    while len(result) != amount:
        totalCount += 1
        result.add(randrange(inclusiveMin, exclusiveMax))
    return result

for i in range(1000):
    getDistinctRandomNumbers(0,1000,632)
print(totalCount*.000001)


# In[23]:

#Possible use as a coefficient of randomness
from random import randrange
totalCount = 0
def biasedRandRange(inclusiveMin, exclusiveMax):
    return int(sqrt(randrange(inclusiveMin, exclusiveMax)*randrange(inclusiveMin, exclusiveMax)))
def getDistinctRandomNumbers(inclusiveMin, exclusiveMax, amount):
    global totalCount
    result = set()
    while len(result) != amount:
        totalCount += 1
        result.add(biasedRandRange(inclusiveMin, exclusiveMax))
    return result

for i in range(100):
    getDistinctRandomNumbers(0,100,62)
print(totalCount*.0001)


# In[24]:

#An implementation of the getRandomElements algorithm
from random import randrange
def shuffle(list):
    n = len(list)
    for i in range(n - 1):
        j = randrange(i, n-1) # i <= j < n
        value = list[i]
        list[i] = list[j]
        list[j] = value
def getDistinctRandomNumbers(inclusiveMin, exclusiveMax, amount):
    result = set()
    while len(result) != amount:
        result.add(randrange(inclusiveMin, exclusiveMax))
    return result
def getRandomElements(list, amount):
    length = len(list)
    amount = min(length, amount)
    if amount < 0.63212055882 * length:
        randomNumbers = getDistinctRandomNumbers(0, length, amount)
        result = [None]*amount
        i = 0
        for num in randomNumbers:
            result[i] = list[num]
            i+=1
        return result
    else:
        shuffle(list)
        return list[0:amount]


# In[25]:

#Example of how to use getRandomElements
import string
print(getRandomElements(list(string.ascii_lowercase),5))

