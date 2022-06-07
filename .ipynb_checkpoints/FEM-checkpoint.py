{
 "cells": [
  {
   "cell_type": "code",
   "execution_count": 565,
   "id": "ea858970",
   "metadata": {},
   "outputs": [],
   "source": [
    "%matplotlib inline\n",
    "\n",
    "import numpy as np\n",
    "import matplotlib.pyplot as plt\n",
    "import scipy\n",
    "import scipy.integrate as integrate\n",
    "import numpy.linalg as linalg"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 566,
   "id": "4d0bb69d",
   "metadata": {},
   "outputs": [],
   "source": [
    "#must be greater than 0\n",
    "startR=0\n",
    "endR=5\n",
    "N=10000\n",
    "sites=N+1\n",
    "deltaX=np.abs(startR-endR)/N\n",
    "eps=1\n",
    "n=1\n",
    "a=1"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 567,
   "id": "eaee86b4",
   "metadata": {},
   "outputs": [],
   "source": [
    "def electronDensity(r, a):\n",
    "    return ((1/(np.pi*a**3))*(np.exp(-2*r/a))) if r>0 else ((1/(np.pi*a**3))*(np.exp(2*r/a)))\n",
    "    \n",
    "def linearFunction(r, n):\n",
    "    return 1/(deltaX) *(r-n)+1 if r<n else -1/(deltaX) *(r-n)+1\n",
    "    \n",
    "def integrand(r, eps, n):\n",
    "    return (1/eps)*electronDensity(r,a)*linearFunction(r, n)*r**2\n",
    "\n",
    "def getIntegral(a_1, a_2, eps, n):\n",
    "    return integrate.quad(integrand, a_1, a_2, args=(eps, n))[0]\n",
    "\n",
    "def v(r):\n",
    "    return (1/(4*np.pi*eps*a**2))*(a**2/r-(a**2/r-a)*np.exp(-2*r/a)) if r>0 else (1/(4*np.pi*eps))*(-a**2/r-(-a**2/r-a)*np.exp(2*r/a))"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 447,
   "id": "559eda95",
   "metadata": {},
   "outputs": [],
   "source": [
    "#uniform\n",
    "def constructResultVector():\n",
    "    theList=[getIntegral(startR+(i-1)*deltaX, startR+(i+1)*deltaX, eps, startR+i*deltaX) for i in range(1,sites-1)];\n",
    "    return np.array(theList);\n",
    "\n",
    "def mk(i):\n",
    "    return ((startR+deltaX*i)**3-(startR+deltaX*(i-1))**3)/(3*deltaX**2)\n",
    "\n",
    "def boundary_condition(r):\n",
    "    if(0<=r and r<endR/2):\n",
    "        return 3/(4*np.pi*a*eps)\n",
    "    return v(r)\n",
    "\n",
    "def addBoundaryConditions(b):\n",
    "    b[0]+=boundary_condition(startR)*mk(1)\n",
    "    b[b.size-1]+=boundary_condition(endR)*mk(N)\n",
    "    return b"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 674,
   "id": "b517fb0b",
   "metadata": {},
   "outputs": [],
   "source": [
    "#func gives grid size at each point\n",
    "def generateUniformGrid(N):\n",
    "    delX=(endR-startR)/N\n",
    "    return [(startR+i*delX,startR+(i+1)*delX) for i in range(N)]\n",
    "\n",
    "def riemannEst(a_1, a_2, func):\n",
    "    return func((a_2-a_1)/2+a_1,a)*(a_2-a_1)\n",
    "\n",
    "def adaptivelySubdivide(l, func, eps,n,maxIt):\n",
    "    absError=np.abs(riemannEst(l[0], l[1], func)-(riemannEst(l[0], (l[1]-l[0])/2+l[0], func)+riemannEst((l[1]-l[0])/2+l[0], l[1],func)))\n",
    "    relError=absError/(riemannEst(l[0], l[1], func))\n",
    "    #print(relError)\n",
    "    if(relError<eps or n>maxIt):\n",
    "        return np.array([l])\n",
    "    else:\n",
    "        left=(l[0], (l[1]-l[0])/2+l[0])\n",
    "        right=((l[1]-l[0])/2+l[0], l[1])\n",
    "        temp=(np.concatenate((adaptivelySubdivide(left, func, eps, n+1,maxIt),adaptivelySubdivide(right, func, eps,n+1,maxIt))))\n",
    "        return temp\n",
    "    \n",
    "#adaptively subdivided to func\n",
    "def generateListOfAdaptiveGridCoords(N, func, eps, maxIt):\n",
    "    uni=generateUniformGrid(N)\n",
    "    newList=[];\n",
    "    for i in range(len(uni)):\n",
    "        #adaptively subdivide returns a list\n",
    "        temp=adaptivelySubdivide(uni[i], func, eps, 0, maxIt);\n",
    "        for x in temp:\n",
    "            newList.append(x)\n",
    "    return np.array(newList,dtype=object)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 677,
   "id": "b7ae8078",
   "metadata": {},
   "outputs": [],
   "source": [
    "temp=generateListOfAdaptiveGridCoords(1000, electronDensity, 1e-4, 5)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 678,
   "id": "88e132fe",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "1000"
      ]
     },
     "execution_count": 678,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "len(temp)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 673,
   "id": "d688f4c4",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "array([[0.0, 0.015625],\n",
       "       [0.015625, 0.03125],\n",
       "       [0.03125, 0.046875],\n",
       "       [0.046875, 0.0625],\n",
       "       [0.0625, 0.078125],\n",
       "       [0.078125, 0.09375],\n",
       "       [0.09375, 0.109375],\n",
       "       [0.109375, 0.125],\n",
       "       [0.125, 0.140625],\n",
       "       [0.140625, 0.15625],\n",
       "       [0.15625, 0.171875],\n",
       "       [0.171875, 0.1875],\n",
       "       [0.1875, 0.203125],\n",
       "       [0.203125, 0.21875],\n",
       "       [0.21875, 0.234375],\n",
       "       [0.234375, 0.25],\n",
       "       [0.25, 0.265625],\n",
       "       [0.265625, 0.28125],\n",
       "       [0.28125, 0.296875],\n",
       "       [0.296875, 0.3125],\n",
       "       [0.3125, 0.328125],\n",
       "       [0.328125, 0.34375],\n",
       "       [0.34375, 0.359375],\n",
       "       [0.359375, 0.375],\n",
       "       [0.375, 0.390625],\n",
       "       [0.390625, 0.40625],\n",
       "       [0.40625, 0.421875],\n",
       "       [0.421875, 0.4375],\n",
       "       [0.4375, 0.453125],\n",
       "       [0.453125, 0.46875],\n",
       "       [0.46875, 0.484375],\n",
       "       [0.484375, 0.5],\n",
       "       [0.5, 0.515625],\n",
       "       [0.515625, 0.53125],\n",
       "       [0.53125, 0.546875],\n",
       "       [0.546875, 0.5625],\n",
       "       [0.5625, 0.578125],\n",
       "       [0.578125, 0.59375],\n",
       "       [0.59375, 0.609375],\n",
       "       [0.609375, 0.625],\n",
       "       [0.625, 0.640625],\n",
       "       [0.640625, 0.65625],\n",
       "       [0.65625, 0.671875],\n",
       "       [0.671875, 0.6875],\n",
       "       [0.6875, 0.703125],\n",
       "       [0.703125, 0.71875],\n",
       "       [0.71875, 0.734375],\n",
       "       [0.734375, 0.75],\n",
       "       [0.75, 0.765625],\n",
       "       [0.765625, 0.78125],\n",
       "       [0.78125, 0.796875],\n",
       "       [0.796875, 0.8125],\n",
       "       [0.8125, 0.828125],\n",
       "       [0.828125, 0.84375],\n",
       "       [0.84375, 0.859375],\n",
       "       [0.859375, 0.875],\n",
       "       [0.875, 0.890625],\n",
       "       [0.890625, 0.90625],\n",
       "       [0.90625, 0.921875],\n",
       "       [0.921875, 0.9375],\n",
       "       [0.9375, 0.953125],\n",
       "       [0.953125, 0.96875],\n",
       "       [0.96875, 0.984375],\n",
       "       [0.984375, 1.0],\n",
       "       [1.0, 1.015625],\n",
       "       [1.015625, 1.03125],\n",
       "       [1.03125, 1.046875],\n",
       "       [1.046875, 1.0625],\n",
       "       [1.0625, 1.078125],\n",
       "       [1.078125, 1.09375],\n",
       "       [1.09375, 1.109375],\n",
       "       [1.109375, 1.125],\n",
       "       [1.125, 1.140625],\n",
       "       [1.140625, 1.15625],\n",
       "       [1.15625, 1.171875],\n",
       "       [1.171875, 1.1875],\n",
       "       [1.1875, 1.203125],\n",
       "       [1.203125, 1.21875],\n",
       "       [1.21875, 1.234375],\n",
       "       [1.234375, 1.25],\n",
       "       [1.25, 1.265625],\n",
       "       [1.265625, 1.28125],\n",
       "       [1.28125, 1.296875],\n",
       "       [1.296875, 1.3125],\n",
       "       [1.3125, 1.328125],\n",
       "       [1.328125, 1.34375],\n",
       "       [1.34375, 1.359375],\n",
       "       [1.359375, 1.375],\n",
       "       [1.375, 1.390625],\n",
       "       [1.390625, 1.40625],\n",
       "       [1.40625, 1.421875],\n",
       "       [1.421875, 1.4375],\n",
       "       [1.4375, 1.453125],\n",
       "       [1.453125, 1.46875],\n",
       "       [1.46875, 1.484375],\n",
       "       [1.484375, 1.5],\n",
       "       [1.5, 1.515625],\n",
       "       [1.515625, 1.53125],\n",
       "       [1.53125, 1.546875],\n",
       "       [1.546875, 1.5625],\n",
       "       [1.5625, 1.578125],\n",
       "       [1.578125, 1.59375],\n",
       "       [1.59375, 1.609375],\n",
       "       [1.609375, 1.625],\n",
       "       [1.625, 1.640625],\n",
       "       [1.640625, 1.65625],\n",
       "       [1.65625, 1.671875],\n",
       "       [1.671875, 1.6875],\n",
       "       [1.6875, 1.703125],\n",
       "       [1.703125, 1.71875],\n",
       "       [1.71875, 1.734375],\n",
       "       [1.734375, 1.75],\n",
       "       [1.75, 1.765625],\n",
       "       [1.765625, 1.78125],\n",
       "       [1.78125, 1.796875],\n",
       "       [1.796875, 1.8125],\n",
       "       [1.8125, 1.828125],\n",
       "       [1.828125, 1.84375],\n",
       "       [1.84375, 1.859375],\n",
       "       [1.859375, 1.875],\n",
       "       [1.875, 1.890625],\n",
       "       [1.890625, 1.90625],\n",
       "       [1.90625, 1.921875],\n",
       "       [1.921875, 1.9375],\n",
       "       [1.9375, 1.953125],\n",
       "       [1.953125, 1.96875],\n",
       "       [1.96875, 1.984375],\n",
       "       [1.984375, 2.0],\n",
       "       [2.0, 2.015625],\n",
       "       [2.015625, 2.03125],\n",
       "       [2.03125, 2.046875],\n",
       "       [2.046875, 2.0625],\n",
       "       [2.0625, 2.078125],\n",
       "       [2.078125, 2.09375],\n",
       "       [2.09375, 2.109375],\n",
       "       [2.109375, 2.125],\n",
       "       [2.125, 2.140625],\n",
       "       [2.140625, 2.15625],\n",
       "       [2.15625, 2.171875],\n",
       "       [2.171875, 2.1875],\n",
       "       [2.1875, 2.203125],\n",
       "       [2.203125, 2.21875],\n",
       "       [2.21875, 2.234375],\n",
       "       [2.234375, 2.25],\n",
       "       [2.25, 2.265625],\n",
       "       [2.265625, 2.28125],\n",
       "       [2.28125, 2.296875],\n",
       "       [2.296875, 2.3125],\n",
       "       [2.3125, 2.328125],\n",
       "       [2.328125, 2.34375],\n",
       "       [2.34375, 2.359375],\n",
       "       [2.359375, 2.375],\n",
       "       [2.375, 2.390625],\n",
       "       [2.390625, 2.40625],\n",
       "       [2.40625, 2.421875],\n",
       "       [2.421875, 2.4375],\n",
       "       [2.4375, 2.453125],\n",
       "       [2.453125, 2.46875],\n",
       "       [2.46875, 2.484375],\n",
       "       [2.484375, 2.5],\n",
       "       [2.5, 2.515625],\n",
       "       [2.515625, 2.53125],\n",
       "       [2.53125, 2.546875],\n",
       "       [2.546875, 2.5625],\n",
       "       [2.5625, 2.578125],\n",
       "       [2.578125, 2.59375],\n",
       "       [2.59375, 2.609375],\n",
       "       [2.609375, 2.625],\n",
       "       [2.625, 2.640625],\n",
       "       [2.640625, 2.65625],\n",
       "       [2.65625, 2.671875],\n",
       "       [2.671875, 2.6875],\n",
       "       [2.6875, 2.703125],\n",
       "       [2.703125, 2.71875],\n",
       "       [2.71875, 2.734375],\n",
       "       [2.734375, 2.75],\n",
       "       [2.75, 2.765625],\n",
       "       [2.765625, 2.78125],\n",
       "       [2.78125, 2.796875],\n",
       "       [2.796875, 2.8125],\n",
       "       [2.8125, 2.828125],\n",
       "       [2.828125, 2.84375],\n",
       "       [2.84375, 2.859375],\n",
       "       [2.859375, 2.875],\n",
       "       [2.875, 2.890625],\n",
       "       [2.890625, 2.90625],\n",
       "       [2.90625, 2.921875],\n",
       "       [2.921875, 2.9375],\n",
       "       [2.9375, 2.953125],\n",
       "       [2.953125, 2.96875],\n",
       "       [2.96875, 2.984375],\n",
       "       [2.984375, 3.0],\n",
       "       [3.0, 3.015625],\n",
       "       [3.015625, 3.03125],\n",
       "       [3.03125, 3.046875],\n",
       "       [3.046875, 3.0625],\n",
       "       [3.0625, 3.078125],\n",
       "       [3.078125, 3.09375],\n",
       "       [3.09375, 3.109375],\n",
       "       [3.109375, 3.125],\n",
       "       [3.125, 3.140625],\n",
       "       [3.140625, 3.15625],\n",
       "       [3.15625, 3.171875],\n",
       "       [3.171875, 3.1875],\n",
       "       [3.1875, 3.203125],\n",
       "       [3.203125, 3.21875],\n",
       "       [3.21875, 3.234375],\n",
       "       [3.234375, 3.25],\n",
       "       [3.25, 3.265625],\n",
       "       [3.265625, 3.28125],\n",
       "       [3.28125, 3.296875],\n",
       "       [3.296875, 3.3125],\n",
       "       [3.3125, 3.328125],\n",
       "       [3.328125, 3.34375],\n",
       "       [3.34375, 3.359375],\n",
       "       [3.359375, 3.375],\n",
       "       [3.375, 3.390625],\n",
       "       [3.390625, 3.40625],\n",
       "       [3.40625, 3.421875],\n",
       "       [3.421875, 3.4375],\n",
       "       [3.4375, 3.453125],\n",
       "       [3.453125, 3.46875],\n",
       "       [3.46875, 3.484375],\n",
       "       [3.484375, 3.5],\n",
       "       [3.5, 3.515625],\n",
       "       [3.515625, 3.53125],\n",
       "       [3.53125, 3.546875],\n",
       "       [3.546875, 3.5625],\n",
       "       [3.5625, 3.578125],\n",
       "       [3.578125, 3.59375],\n",
       "       [3.59375, 3.609375],\n",
       "       [3.609375, 3.625],\n",
       "       [3.625, 3.640625],\n",
       "       [3.640625, 3.65625],\n",
       "       [3.65625, 3.671875],\n",
       "       [3.671875, 3.6875],\n",
       "       [3.6875, 3.703125],\n",
       "       [3.703125, 3.71875],\n",
       "       [3.71875, 3.734375],\n",
       "       [3.734375, 3.75],\n",
       "       [3.75, 3.765625],\n",
       "       [3.765625, 3.78125],\n",
       "       [3.78125, 3.796875],\n",
       "       [3.796875, 3.8125],\n",
       "       [3.8125, 3.828125],\n",
       "       [3.828125, 3.84375],\n",
       "       [3.84375, 3.859375],\n",
       "       [3.859375, 3.875],\n",
       "       [3.875, 3.890625],\n",
       "       [3.890625, 3.90625],\n",
       "       [3.90625, 3.921875],\n",
       "       [3.921875, 3.9375],\n",
       "       [3.9375, 3.953125],\n",
       "       [3.953125, 3.96875],\n",
       "       [3.96875, 3.984375],\n",
       "       [3.984375, 4.0],\n",
       "       [4.0, 4.015625],\n",
       "       [4.015625, 4.03125],\n",
       "       [4.03125, 4.046875],\n",
       "       [4.046875, 4.0625],\n",
       "       [4.0625, 4.078125],\n",
       "       [4.078125, 4.09375],\n",
       "       [4.09375, 4.109375],\n",
       "       [4.109375, 4.125],\n",
       "       [4.125, 4.140625],\n",
       "       [4.140625, 4.15625],\n",
       "       [4.15625, 4.171875],\n",
       "       [4.171875, 4.1875],\n",
       "       [4.1875, 4.203125],\n",
       "       [4.203125, 4.21875],\n",
       "       [4.21875, 4.234375],\n",
       "       [4.234375, 4.25],\n",
       "       [4.25, 4.265625],\n",
       "       [4.265625, 4.28125],\n",
       "       [4.28125, 4.296875],\n",
       "       [4.296875, 4.3125],\n",
       "       [4.3125, 4.328125],\n",
       "       [4.328125, 4.34375],\n",
       "       [4.34375, 4.359375],\n",
       "       [4.359375, 4.375],\n",
       "       [4.375, 4.390625],\n",
       "       [4.390625, 4.40625],\n",
       "       [4.40625, 4.421875],\n",
       "       [4.421875, 4.4375],\n",
       "       [4.4375, 4.453125],\n",
       "       [4.453125, 4.46875],\n",
       "       [4.46875, 4.484375],\n",
       "       [4.484375, 4.5],\n",
       "       [4.5, 4.515625],\n",
       "       [4.515625, 4.53125],\n",
       "       [4.53125, 4.546875],\n",
       "       [4.546875, 4.5625],\n",
       "       [4.5625, 4.578125],\n",
       "       [4.578125, 4.59375],\n",
       "       [4.59375, 4.609375],\n",
       "       [4.609375, 4.625],\n",
       "       [4.625, 4.640625],\n",
       "       [4.640625, 4.65625],\n",
       "       [4.65625, 4.671875],\n",
       "       [4.671875, 4.6875],\n",
       "       [4.6875, 4.703125],\n",
       "       [4.703125, 4.71875],\n",
       "       [4.71875, 4.734375],\n",
       "       [4.734375, 4.75],\n",
       "       [4.75, 4.765625],\n",
       "       [4.765625, 4.78125],\n",
       "       [4.78125, 4.796875],\n",
       "       [4.796875, 4.8125],\n",
       "       [4.8125, 4.828125],\n",
       "       [4.828125, 4.84375],\n",
       "       [4.84375, 4.859375],\n",
       "       [4.859375, 4.875],\n",
       "       [4.875, 4.890625],\n",
       "       [4.890625, 4.90625],\n",
       "       [4.90625, 4.921875],\n",
       "       [4.921875, 4.9375],\n",
       "       [4.9375, 4.953125],\n",
       "       [4.953125, 4.96875],\n",
       "       [4.96875, 4.984375],\n",
       "       [4.984375, 5.0]], dtype=object)"
      ]
     },
     "execution_count": 673,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "temp"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 449,
   "id": "f30fbfcc",
   "metadata": {},
   "outputs": [],
   "source": [
    "def constructMatrix():\n",
    "    temp=np.zeros((sites, sites));\n",
    "    for i in range(sites-1):\n",
    "        if(i>0):\n",
    "            temp[i, i-1]=-mk(i)\n",
    "        if(i>=0):\n",
    "            temp[i, i]=mk(i+1)+mk(i)\n",
    "        if(i+1<sites):\n",
    "            temp[i, i+1]=-mk(i+1)\n",
    "    temp=np.delete(temp, [0, sites-1], axis=1)\n",
    "    temp=np.delete(temp, [0, sites-1], axis=0)\n",
    "    return temp"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 450,
   "id": "85e98221",
   "metadata": {},
   "outputs": [],
   "source": [
    "b=constructResultVector()\n",
    "b=addBoundaryConditions(b)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 451,
   "id": "eb6e4d33",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "(9999,)"
      ]
     },
     "execution_count": 451,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "np.shape(b)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 452,
   "id": "c8194ded",
   "metadata": {},
   "outputs": [],
   "source": [
    "M=constructMatrix()"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 453,
   "id": "bfb46232",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "(9999, 9999)"
      ]
     },
     "execution_count": 453,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "np.shape(M)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 454,
   "id": "289ccf60",
   "metadata": {},
   "outputs": [],
   "source": [
    "sols=linalg.solve(M, b)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 466,
   "id": "ae28c95c",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "x1=[startR+i*deltaX for i in range(b.size)]\n",
    "plt.scatter(x1,sols,c=\"red\")\n",
    "x = np.linspace(startR,endR,1000)\n",
    "y=[v(r) for r in x]\n",
    "#plt.scatter(x, y)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 467,
   "id": "722ccb3f",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "<matplotlib.collections.PathCollection at 0x1243ed760>"
      ]
     },
     "execution_count": 467,
     "metadata": {},
     "output_type": "execute_result"
    },
    {
     "data": {
      "image/png": "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\n",
      "text/plain": [
       "<Figure size 432x288 with 1 Axes>"
      ]
     },
     "metadata": {
      "needs_background": "light"
     },
     "output_type": "display_data"
    }
   ],
   "source": [
    "plt.scatter(x, y)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 456,
   "id": "f76d3d46",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "0.238732414637843"
      ]
     },
     "execution_count": 456,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "3/(4*np.pi*a*eps)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 457,
   "id": "8865dcce",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "array([0.01591851, 0.01591864, 0.01591876, ..., 0.01591876, 0.01591864,\n",
       "       0.01591851])"
      ]
     },
     "execution_count": 457,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "sols"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 458,
   "id": "d57f3b93",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "0.015918384558484624"
      ]
     },
     "execution_count": 458,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "y[np.shape(y)[0]-1]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 459,
   "id": "e41fc932",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "1000"
      ]
     },
     "execution_count": 459,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "np.shape(y)[0]"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 460,
   "id": "e9b8ce7f",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "0.015918384558484624"
      ]
     },
     "execution_count": 460,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "boundary_condition(startR)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 461,
   "id": "6cb42193",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "0.015918384558484624"
      ]
     },
     "execution_count": 461,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "boundary_condition(endR)"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": 462,
   "id": "f27ff7bb",
   "metadata": {},
   "outputs": [
    {
     "data": {
      "text/plain": [
       "2"
      ]
     },
     "execution_count": 462,
     "metadata": {},
     "output_type": "execute_result"
    }
   ],
   "source": [
    "lol=np.array([0,1])\n",
    "lol.size"
   ]
  },
  {
   "cell_type": "code",
   "execution_count": null,
   "id": "2564732f",
   "metadata": {},
   "outputs": [],
   "source": []
  }
 ],
 "metadata": {
  "kernelspec": {
   "display_name": "Python 3 (ipykernel)",
   "language": "python",
   "name": "python3"
  },
  "language_info": {
   "codemirror_mode": {
    "name": "ipython",
    "version": 3
   },
   "file_extension": ".py",
   "mimetype": "text/x-python",
   "name": "python",
   "nbconvert_exporter": "python",
   "pygments_lexer": "ipython3",
   "version": "3.9.12"
  }
 },
 "nbformat": 4,
 "nbformat_minor": 5
}
