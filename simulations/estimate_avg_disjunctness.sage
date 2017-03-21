import time
load rs.sage
#################################################################
### This file estimates the average disjunctness of:
### (1) Our matrix (first cosets)
### (2) Random columns of an RS matrix
### (3) Random matrix
##################################################################

## N is the number of specimens, m is the number of layers to use.  In some cases, the basic matrix will already have been computed, you can pass this in in basemat.
def make_our_rs_mat(N=480,m=6,q=16,basemat=None):
	if basemat is None:
		deg = int(log( N, q ))  # make sure we have enough columns to cover N
		# this is the ordering that will order polynomials like "0, 1, 2, .. ,15, x, x+1, x+2, .. ,x+15, 2x, 2x + 1, ...."
		C1,image,polys = create_rs_code( q, deg, concatenated=True, orderoption="lex", orderfieldelts="lex")
	else:
		C1 = basemat
	return C1[ :m*q, :N ] # return the first N columns (in our coset order) and the first 6 layers.  

def make_rand_rs_mat(N=480,m=6,q=16,basemat=None):
	if basemat is None:
		deg = int(log( N, q )) 
		C1,image,polys = create_rs_code( q, deg, concatenated=True, orderoption="lex", orderfieldelts="lex")
	else:
		C1 = basemat
	# make a random subset 
	choices = range(C1.ncols())
	shuffle(choices)
	return C1[ :m*q, choices[:N] ] # return a random N columns and the first m layers.  
	
def make_rand_mat(N=480, m=6,q=16):
	ret = matrix( [ [0 for i in range(N) ] for j in range(m*q) ] )
	for i in range(N):
		T = random_subset( m, range(q*m) ) # a random subset of 1's -- there should be m of them in each column
		for j in T:
			ret[j,i] = 1
	return ret

def random_subset(k, S):
	T = S[:]
	shuffle(T)
	return T[:k]


######  ESTIMATE DISJUNCTNESS ########################

## true if the column i and the set S is good for M
def isGoodCol( i, S, M ):
	if i in S:
		print "BAD!"
		return None
	for l in range(M.nrows()):
		if M[l,i] == 1:
			covered = False
			for j in S:
				if M[l,j] == 1:
					covered = True
					break
			if not covered:
				return True
	return False
			

## true if the set S of columns is good for M
def isGood( S, M ):
	for i in range(M.ncols()):
		if i in S:
			continue
		if not isGoodCol( i, S, M ):
			return False
	return True
		
# returns the probability of success for k trials
def est_avg_disjunctness( M, k, trials=5000 ):
	good_count = 0
	for t in range(trials):
		S = random_subset(k, range(M.ncols()))
		if isGood(S, M):
			good_count += 1
	p = good_count / trials
	stderr = getstderr(p, trials)
	return float(p), float(stderr)

# when we have a bunch of 0/1 random variables, the standard error is determined by the sample mean 
# and the number of trials
def getstderr( p , trials):
	stderr = sqrt( p* (1 - p)^2 + (1 -p)* p^2 ) / sqrt(trials)
	return float(stderr)


### output our data

def do_some_trials(M, trials, ks ):
	ests = {}
	errs = {}
	flag = False
	for k in ks:
		if flag:
			est = 0
			err = 0
		else:
			est,err = est_avg_disjunctness( M, k, trials=trials )
		ests[k] = est
		errs[k] = err
		if est == 0:
			flag = True  # once I get to zero, stop wasting time
	return ests, errs
		


def generate_tables(trials=1000, Ns=[96*4,96*5,96*6], fbasename="avg_disj",q=16,randtrials=100,ks=range(20),dorand=True):
	uniq = str( time.asctime() )
	# in most cases, the same matrix will do for all the constructions.  In that case, don't do it lots of times.
	if all( [ int( log(N,16) ) == int( log(Ns[0],16) ) for N in Ns ] ):
		basemat, image,polys  = create_rs_code( q, int(log(Ns[0],16)) , concatenated=True, orderoption="lex", orderfieldelts="lex")
	else:
		basemat = None
	for N in Ns:
		print "N=",N
		Mours = make_our_rs_mat(N=N, m=6,q=16,basemat=basemat)
		## do it once for our deterministic matrix
		ests, errs = do_some_trials(Mours, trials, ks)
		
		if dorand:
			## now do it a bunch of times for random matrices
			eststmp1 = {}
			eststmp2 = {}
			for s in range(randtrials):
				Mrand1 = make_rand_rs_mat(N=N,m=6,q=16,basemat=basemat)
				estsrand1, errsrand1 = do_some_trials( Mrand1, trials, ks)
				for k in estsrand1.keys():
					if k not in eststmp1:
						eststmp1[k] = []
					eststmp1[k].append( estsrand1[k] )
				Mrand2 = make_rand_mat(N=N,m=6,q=16)
				estsrand2, errsrand2 = do_some_trials( Mrand2, trials, ks )
				for k in estsrand2.keys():
					if k not in eststmp2:
						eststmp2[k] = []
					eststmp2[k].append( estsrand2[k] )
	
			# figure out standard errors and stuff from the random trials
			ests1 = {}
			ests2 = {}
			errs1 = {}
			errs2 = {}
			for k in eststmp1.keys():
				avg = sum( eststmp1[k] ) / len( eststmp1[k] )
				ests1[k] = float(avg)
				err = sqrt(sum([  (x - avg)^2 for x in eststmp1[k] ]))/sqrt(len(eststmp1[k]))
				errs1[k] = float(err)
			for k in eststmp2.keys():
				avg = sum( eststmp2[k] ) / len( eststmp2[k] )
				ests2[k] = float(avg)
				err = sqrt(sum([  (x - avg)^2 for x in eststmp2[k] ]))/sqrt(len(eststmp2[k]))
				errs2[k] = float(err)
		
		# finally, write everything to a file
		F = open( fbasename + "_" + uniq  + "_" + str(N) + ".dat" , "w")
		for k in ks:
			F.writelines(str(k) + "\t")
			F.writelines( str(ests[k]) + "\t" + str(errs[k]) + "\t"  )  # our matrix (error bars are redundant)
			if dorand:
				F.writelines( str(ests1[k]) + "\t" + str(errs1[k]) + "\t" + str( getstderr( ests1[k] , trials) ) + "\t") # rand1
				F.writelines( str(ests2[k]) + "\t" + str(errs2[k]) + "\t" + str( getstderr( ests2[k] , trials) ) + "\t") # rand2
			F.writelines("\n")
		F.close()

