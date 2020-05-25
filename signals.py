"""
====================================================================================
signals
====================================================================================
Description: signals available on wavelab
====================================================================================
"""
import sys
import numpy as np

def exception_handler(exception_type, exception, traceback):
	print "%s: %s" % (exception_type.__name__, exception)
sys.excepthook = exception_handler


def generateHeaviSine(N):
	t = np.linspace(1, N, N) / float(N)
	y = 4.0*np.sin(4.0*np.pi*t) - np.sign(t-0.3) - np.sign(0.72-t)
	return t, y

def generateBlocks(N):
	pos = [0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81]
	hgt = [4., -5., 3., -4., 5., -4.2, 2.1, 4.3, -3.1, 2.1, -4.2]
	t = np.linspace(1, N, N) / float(N)
	y = np.zeros(N)
	for p, h in zip(pos, hgt):
		y += (h/2.0) * (1.0 + np.sign(t - p))
	return t, y

def generateBumps(N):
	pos = [0.10, 0.13, 0.15, 0.23, 0.25, 0.40, 0.44, 0.65, 0.76, 0.78, 0.81]
	hgt = [4., 5., 3., 4., 5., 4.2, 2.1, 4.3, 3.1, 5.1, 4.2]
	wth = [0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005]
	t = np.linspace(1, N, N) / float(N)
	y = np.zeros(N)
	for p, h, w in zip(pos, hgt, wth):
		y += h * ((1.0+np.abs((t-p)/w))**-4)
	return t, y

def generateDoppler(N, epsilon = 0.05):
	t = np.linspace(1, N, N) / float(N)
	y = np.sqrt(t *(1.0-t)) * np.sin((2.0*np.pi*(1.0+epsilon))/(t+epsilon))
	return t, y

def generateQuadchirp(N):
	t = np.linspace(1, N, N) / float(N)
	y = np.sin((np.pi/3.0)*N*t*t*t)
	return t, y

def generateMishmash(N):
	# QuadChirp + LinChirp + HiSine
	t  = np.linspace(1, N, N) / float(N)
	y  = np.sin((np.pi/3.0)*N*t*t*t) + np.sin(0.125*np.pi*N*t*t) + np.sin(0.6902*np.pi*N*t)
	return t, y

def generateRamp(N):
	t = np.linspace(1, N, N) / float(N)
	y = t - (t>=0.37)
	return t, y 

def generateCusp(N):
	t = np.linspace(1, N, N) / float(N)
	y = np.sqrt(np.abs(t-0.37))
	return t, y

def generateSing(N):
	t = np.linspace(1, N, N) / float(N)
	k = np.floor(N*0.37)
	y = 1.0 / np.abs(t-(k+0.5)/N)
	return t, y

def generateHiSine(N):
	t = np.linspace(1, N, N) / float(N)
	y = np.sin(0.6902*np.pi*N*t)
	return t, y

def generateLoSine(N):
	t = np.linspace(1, N, N) / float(N)
	y = np.sin(0.3333*np.pi*N*t)
	return t, y

def generateLinChirp(N):
	t = np.linspace(1, N, N) / float(N)
	y = np.sin(N*0.500*np.pi*t*t)
	#y = np.sin(N*0.125*np.pi*t*t)
	return t, y

def generateTwoChirp(N):
	t = np.linspace(1, N, N) / float(N)
	y = np.sin(np.pi*N*t*t) + np.sin((np.pi/3.0)*N*t*t)
	return t, y

def generateWernerSorrows(N):
	t = np.linspace(1, N, N) / float(N)
	y  = np.sin(np.pi*(N/2.0)*t*t*t)
	y += np.sin(0.6902*np.pi*N*t)
	y += np.sin(np.pi*N*t*t)
	pos = [0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81]
	hgt = [4., 5., 3., 4., 5., 4.2, 2.1, 4.3, 3.1, 5.1, 4.2]
	wth = [0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005]
	for p, h, w in zip(pos, hgt, wth):
		y += h / (1.0+np.abs((t-p)/w))**4
	return t, y

def generateLeopold(N):
	t = np.linspace(1, N, N) / float(N)
	y = (t == (np.floor(0.37*N)/float(N))).astype(np.float)
	return t, y

def generateRiemann(N):
	t = np.linspace(1, N, N) / float(N)
	sqn = int(np.round(np.sqrt(N)))
	y = t * 0
	for i in xrange(sqn):
		y[((i+1)**2)-1] = 1.0 / (i+1)
	y = np.real(np.fft.ifft(y))
	return t, y

def generateHypChirps(N):
	alpha = 15.0 * N * np.pi / 1024
	beta  = 5.0 * N * np.pi / 1024
	t  	= np.arange(1.001, N+0.001+1, 1) / float(N) 
	f1 = np.sin(alpha / (0.8 - t)) * (0.1 < t) * (t < 0.68)
	f2 = np.sin(beta / (0.8 - t)) * (0.1 < t) * (t < 0.75)
	M  = int(round(0.65 * N))
	P = int(np.floor(M / 4))
	enveloppe = np.ones(M)
	enveloppe[1:P] = (1.0 + np.sin(-np.pi / 2.0 + (np.arange(1,P,1) - np.ones(P-1)) / (P-1) * np.pi)) / 2.0
	enveloppe[M-P+1:M] = enveloppe[P:1:-1]
	env = np.zeros(N)
	env[int(np.ceil(N / 10)):int(M+np.ceil(N / 10)-1)] = enveloppe[1:M]
	y = (f1+f2) * env
	return t, y

def generateLinChirps(N):
	b 	= 100.0 * N * np.pi / 1024
	a 	= 250.0 * N * np.pi / 1024
	t 	= np.linspace(1, N, N) / float(N)
	A1 	= np.sqrt((t - 1.0 / N) * (1.0 - t))
	y	= A1 * (np.cos(a * t**2) + np.cos(b * t + a * t**2))
	return t, y

def generateChirps(N):
	t 	= (np.linspace(1, N, N) / float(N)) * 10.0 * np.pi
	f1 	= np.cos(t*t*N/1024.0)
	a 	= 30.0 * N / 1024.0
	t 	= (np.linspace(1, N, N) / float(N)) * np.pi  
	f2 	= np.cos(a*(t*t*t))
	f2 	= [f2[i] for i in xrange(N-1, -1,-1)]
	ix 	= (np.arange(-N, N+1, 1) / float(N)) * 20.0
	g 	= np.exp(-ix**2 * 4.0 * N / 1024.0)
	i1 	= np.arange(N/2 + 1, N/2 + N + 1, 1)
	i2 	= np.arange(N/8 + 1, N / 8 + N + 1, 1)
	j  	= np.linspace(1, N, N) / float(N)
	f3 	= g[i1-1] * np.cos(50.0*np.pi*j*N/1024.0)
	f4 	= g[i2-1] * np.cos(350.0*np.pi*j*N/1024.0)
	y 	= f1 + f2 + f3 + f4
	enveloppe = np.ones(N)
	enveloppe[0:N/8] = (1.0 + np.sin(-np.pi / 2.0 + (np.arange(1, N/8+1, 1) - np.ones(N/8)) / (N/8 - 1) * np.pi))/2.0
	enveloppe[7*N/8:N] = [enveloppe[i] for i in xrange(N/8-1, -1, -1)]
	y 	= y * enveloppe
	return t, y

def generatePieceRegular(N):
	sig1 = -15.0 * generateBumps(N)[1]
	t = np.arange(1, np.fix(N/12)+1, 1) / np.fix(N/12)
	sig2 = -np.exp(4.0 * t)
	t = np.arange(1, np.fix(N/7)+1, 1) / np.fix(N/7)
	sig5 = np.exp(4.0 * t) - np.exp(4.0)
	t = np.arange(1, np.fix(N/3)+1, 1) / np.fix(N/3)
	sigma = 6.0 / 40.0
	sig6 = -70.0 * np.exp(-((t-1.0/2.0) * (t-1.0/2.0)) / (2.0 * sigma**2))
	sig = np.zeros(N)
	sig[0:int(np.fix(N/7))] = sig6[0:int(np.fix(N/7))]
	sig[int(np.fix(N/7)):int(np.fix(N/5))] = 0.5 * sig6[int(np.fix(N/7)):int(np.fix(N/5))]
	sig[int(np.fix(N/5)):int(np.fix(N/3))] = sig6[int(np.fix(N/5)):int(np.fix(N/3))]
	sig[int(np.fix(N/3)):int(np.fix(N/2))] = sig1[int(np.fix(N/3)):int(np.fix(N/2))]
	sig[int(np.fix(N/2)):int(np.fix(N/2)+np.fix(N/12))] = sig2
	sig[int(np.fix(N/2)+2*np.fix(N/12)-1):int(np.fix(N/2)+np.fix(N/12)-1):-1] = sig2
	sig[int(np.fix(N/2)+2*np.fix(N/12)+np.fix(N/20)):int(np.fix(N/2)+2*np.fix(N/12)+3*np.fix(N/20))] = -np.ones(int(np.fix(N/2)+2*np.fix(N/12)+3*np.fix(N/20)-np.fix(N/2)-2*np.fix(N/12)-np.fix(N/20))) * 25
	k = np.fix(N/2) + 2.0 * np.fix(N/12) + 3.0 * np.fix(N/20)
	sig[int(k):int(k+np.fix(N/7))] = sig5
	diff = int(N - 5.0 * np.fix(N/5))

	sig[int(5*np.fix(N/5)):N] = [sig[i] for i in xrange(diff-1, -1, -1)]
	bias = np.sum(sig) / float(N)
	y = bias-sig
	t = np.linspace(1, N, N) / float(N)
	return t, y

def generatePiecePolynomial(N):
	t = np.arange(1, np.fix(N/5)+1, 1) / np.fix(N/5)
	sig1 = 20.0 * (t**3 + t**2 + 4)
	sig3 = 40.0 * (2.0 * t**3 + t) + 100
	sig2 = 10.0 * t**3 + 45
	sig4 = 16.0 * t**2 + 8.0*t + 16
	sig5 = 20.0 * (t + 4)
	sig = np.zeros(N)
	sig[0:int(np.fix(N/5))] = sig1
	sig[int(2*np.fix(N/5)-1):int(np.fix(N/5)-1):-1] = sig2
	sig[int(2*np.fix(N/5)):int(3*np.fix(N/5))] = sig3
	sig[int(3*np.fix(N/5)):int(4*np.fix(N/5))] = sig4
	sig[int(4*np.fix(N/5)):int(5*np.fix(N/5))] = np.fliplr([sig5])[0]
	diff = int(N - 5.0 * np.fix(N/5))
	sig[int(5*np.fix(N/5)):N] = [sig[i] for i in xrange(diff-1,-1,-1)]
	sig[int(np.fix(N/20)):int(np.fix(N/20)+np.fix(N/10))] = np.ones(int(np.fix(N/10))) * 10
	sig[int(N-np.fix(N/10)):int(N+np.fix(N/20)-np.fix(N/10))] = np.ones(int(np.fix(N/20))) * 150
	bias = np.sum(sig) / float(N)
	y = sig - bias
	t = np.linspace(1, N, N) / float(N)
	return t, y

def generateGabor():
	# two modulated Gabor functions in Mallat's book
	N = 512
	t = (np.linspace(-N, N, 2*N+1) * 5.0) / float(N)
	j = np.linspace(1, N, N) / float(N)
	g = np.exp(-t**2 * 20.0)
	i1 = [i for i in xrange(2*N/4, 2*N/4+N, 1)]
	i2 = [i for i in xrange(N/4, N/4+N, 1)]
	sig1 = 3.0 * g[i1] * np.exp(1j*N/16.0*np.pi*j)
	sig2 = 3.0 * g[i2] * np.exp(1j*N/4.0*np.pi*j)
	y = sig1 + sig2
	t = np.linspace(1, N, N) / float(N)
	return t, y

def generateSineoneoverx(): 
	# sin(1/x) in Mallat's book
	N = 1024
	i1 = np.array([float(i) for i in xrange(-N+1, N+1, 1)])
	i1[N-1] = 1.0 / 100.0
	i1 = i1 / (N-1)
	sig = np.sin(1.5 / i1)
	y = sig[512:1536]
	t = np.linspace(1, N, N)  / float(N)
	return t, y

def generateCusp2():
	N = 64
	i1 = np.linspace(1, N, N) / float(N)
	x = (1.0-np.sqrt(i1)) + i1/2.0 - 0.5
	M = 8*N
	y = np.zeros(M)
	y[int(M-1.5*N):int(M-0.5*N)] = x
	y[int(M-2.5*N+1):int(M-1.5*N+1)] = [x[i] for i in xrange(len(x)-1, -1, -1)]
	y[int(3*N):int(3*N+N)] = 0.5 * np.ones(N)
	t = np.linspace(1, M, M) / float(M)
	return t, y

def generatePositiveTriangularWave(t0 = -4, tf = 4, dt = 0.01, A = 1.0, T = 2.0):
	t = np.arange(t0, tf+dt, dt)
	y = A*(T*np.abs(np.rint((1.0/T)*t)-((1.0/T)*t)))
	return t, y

def generateSawtoothWave(t0 = -4, tf = 4, dt = 0.01, A = 1.0, T = 2.0, phase = 0, reverse = False):
	t = (np.arange(t0, tf+dt, dt) / T) + phase
	if reverse:
		A = -A
	y = 2*A*(t-np.floor(t))-A
	t = (t*T) - phase
	return t, y

def generateTriangularWave(t0 = -4, tf = 4, dt = 0.01, A = 2.0, T = 2.0, phase = 0):
	t = np.arange(t0, tf+dt, dt) + phase
	y = A*(4.0/T * (t - T/2 * np.floor((2*t)/T + 1.0/2.0))*(-1.0)**np.floor((2*t)/T + 1.0/2.0))
	t = t - phase
	return t, y

def makeSignal(signalName, N = 1024, epsilon = 0.05):
	"""
	Creates artificial test signal identical to the test signals proposed in WaveLab

	Usage
	-----
	t, y = makeSignal(signalName, N)

	Parameters
	----------
	signalName : string
		Name of the desired signal. Supported values:
			* 'HeaviSine'
			* 'Blocks'
			* 'Bumps'
			* 'Doppler'
			* 'QuadChirp'
			* 'MishMash'
			* 'Ramp'
			* 'Cusp'
			* 'Sing'
			* 'HiSine'
			* 'LoSine'
			* 'LinChirp'
			* 'TwoChirp'
			* 'WernerSorrows' (Heisenberg)
			* 'Leopold' (Kronecker)
			* 'Riemann'
			* 'HypChirps'
			* 'LinChirps'
			* 'Chirps'
			* 'Piece-Regular' (Piece-Wise Smooth)
			* 'Piece-Polynomial' (Piece-Wise 3rd degree polynomial)
			* 'Gabor'
			* 'sineoneoverx'
			* 'Cusp2'
	N : integer, optional (default = 1024)
        Desired signal length.
	epsilon: float, optional (default = 0.05)
		Epsilon value for Doppler function

    Returns
    -------
	t : time, length = N
    y : signal, length = N

    References
    ----------
    WaveLab: http://playfair.stanford.edu/~wavelab/

	Comments
	--------
	Functions not implemented: 'SmoothCusp', 'Gaussian'
    """
	if signalName == 'HeaviSine':
		t, y = generateHeaviSine(N)
	elif signalName == 'Blocks':
		t, y = generateBlocks(N)
	elif signalName == 'Bumps':
		t, y = generateBumps(N)
	elif signalName == 'Doppler':
		t, y = generateDoppler(N, epsilon)
	elif signalName == 'QuadChirp':
		t, y = generateQuadchirp(N)
	elif signalName == 'MishMash':
		t, y = generateMishmash(N)
	elif signalName == 'Ramp':
		t, y = generateRamp(N)
	elif signalName == 'Cusp':
		t, y = generateCusp(N)
	elif signalName == 'Sing':
		t, y = generateSing(N)
	elif signalName == 'HiSine':
		t, y = generateHiSine(N)
	elif signalName == 'LoSine':
		t, y = generateLoSine(N)
	elif signalName == 'LinChirp':
		t, y = generateLinChirp(N)
	elif signalName == 'TwoChirp':
		t, y = generateTwoChirp(N)
	elif signalName == 'WernerSorrows':
		t, y = generateWernerSorrows(N)
	elif signalName == 'Leopold':
		t, y = generateLeopold(N)
	elif signalName == 'Riemann':
		t, y = generateRiemann(N)
	elif signalName == 'HypChirps':
		t, y = generateHypChirps(N)
	elif signalName == 'LinChirps':
		t, y = generateLinChirps(N)
	elif signalName == 'Chirps':
		t, y = generateChirps(N)
	elif signalName == 'Piece-Regular':
		t, y = generatePieceRegular(N)
	elif signalName == 'Piece-Polynomial':
		t, y = generatePiecePolynomial(N)
	elif signalName == 'Gabor':
		t, y = generateGabor()
	elif signalName == 'sineoneoverx':
		t, y = generateSineoneoverx()
	elif signalName == 'Cusp2':
		t, y = generateCusp2()
	else:
		raise ValueError("Unrecognized signal name. Allowable names are: 'HeaviSine', 'Blocks', 'Bumps', 'Doppler', 'QuadChirp', 'MishMash', 'Ramp', 'Cusp', 'Sing', 'HiSine', 'LoSine', 'LinChirp', 'TwoChirp', 'WernerSorrows', 'Leopold', 'Riemann', 'HypChirps', 'LinChirps', 'Chirps', 'Piece-Regular', 'Piece-Polynomial', 'Gabor', 'sineoneoverx' and 'Cusp2'")

	return t, y
