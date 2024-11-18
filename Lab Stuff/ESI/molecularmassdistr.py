import numpy as np
import matplotlib.pyplot as plt
from scipy.special import binom

maxIsoC = 28
maxIsoN = 2
maxIsoO = 3
maxIsoCl = 1
nH = 31

probIsoC = 0.01
probIsoN = 0.003
probIsoO = 0.001
probIsoCl = 0.24


mH = 1
mC = 12
mCIso = 13

mN = 14
mNIso = 15

mO = 16
mOIso = 18

mCl = 35
mClIso = 37
checksum =0
probs = []
masses = []


#PT is expected to only hold approx million ions, throw away any probability for which P < 1e-7



mHtot =mH*nH 
for nCIso in range(maxIsoC+1):
    mCtot = nCIso*mCIso + (maxIsoC-nCIso) * mC
    pC = binom(maxIsoC,nCIso) * probIsoC**nCIso * (1-probIsoC)**(maxIsoC-nCIso)
    for nNIso in range(maxIsoN+1):
        mNtot = nNIso*mNIso + (maxIsoN-nNIso) * mN
        pN= binom(maxIsoN,nNIso) * probIsoN**nNIso * (1-probIsoN)**(maxIsoN-nNIso)
        for nOIso in range(maxIsoO+1):
            mOtot = nOIso*mOIso + (maxIsoO-nOIso) * mO
            pO = binom(maxIsoO,nOIso) * probIsoO**nOIso * (1-probIsoO)**(maxIsoO-nOIso)
            for nClIso in range(maxIsoCl+1):
                mCltot = nClIso*mClIso + (maxIsoCl-nClIso) * mCl
                pCl = binom(maxIsoCl,nClIso) * probIsoCl**nClIso * (1-probIsoCl)**(maxIsoCl-nClIso)
                p = pC*pN*pO*pCl
                M = mCtot+mNtot+mCltot+mOtot+mHtot
                probs.append(p)
                masses.append(M)
print(np.sum(probs))

#Add together isotopes that end up with the same masses
tempProb = []
tempmass = []
probs = np.array(probs)
masses = np.array(masses)

#Determine all the unique masses
UniqueMasses = np.unique(masses)
print(UniqueMasses)
mNew = []
pNew = []
#Add duplicates
for m in UniqueMasses:
    indices = np.where(masses == m)
    mNew.append(m)
    pNew.append(np.sum(probs[indices]))

pNew = np.array(pNew)
mNew = np.array(mNew)
probsFiltered = np.array(pNew)[pNew>1e-7]
massesFiltered = np.array(mNew)[pNew>1e-7]
plt.xlabel('Mass / a.m.u.')
plt.ylabel('Probability')
plt.title('Rhodamine 6G mass distribution')

plt.scatter(massesFiltered,probsFiltered)
plt.show()
